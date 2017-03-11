#include <atomic>
#include <mutex>
#include <utility>
#include <memory>

#include "pmemory.hpp"

namespace aux
{
  template <class T>
  struct stack_elem
  {
      stack_elem(const T& e, stack_elem* n)
      : entry(e), nextptr(n)
      {}

      // getters
      T           elem() const { return entry; }
      stack_elem* next() const { return nextptr; }

      // setters
      void set_next(stack_elem* n)
      {
        nextptr = n;
      }

    private:
      T           entry;
      stack_elem* nextptr;
  };
}

namespace lockfree
{
  template <class _Alloc>
  struct ptr_guard
  {
    ptr_guard(_Alloc alloc, typename _Alloc::value_type* elem)
    : node_allocator(alloc), e(elem)
    {}

    ~ptr_guard()
    {
      node_allocator.deallocate(e, 1);
    }

    _Alloc                             node_allocator;
    typename _Alloc::value_type* const e;
  };

  template < class _Tp,
             class _Alloc = lockfree::pub_scan_manager<_Tp, std::allocator>
           >
  struct stack
  {
      typedef aux::stack_elem<_Tp>          stack_elem;

      typedef typename _Alloc::value_type   _Alloc_value_type;
      typedef std::allocator_traits<_Alloc> _OrigAlloc_traits;
      typedef typename _OrigAlloc_traits::template rebind_alloc<_Tp>        _Tp_alloc_type;
      typedef typename _OrigAlloc_traits::template rebind_alloc<stack_elem> _Node_alloc_type;
      typedef std::allocator_traits<_Tp_alloc_type>                         _Alloc_traits;
      typedef typename _Node_alloc_type::pinguard                           PinGuard;
      typedef _Alloc                                                        allocator_type;

      explicit
      stack(const allocator_type& alloc = allocator_type())
      : top(nullptr), nodeAlloc(alloc)
      {}

      void push(const _Tp& e)
      {
        _Node_alloc_type alloc   = get_allocator();
        stack_elem*      newelem = alloc.allocate(1);
        stack_elem*      currtop = nullptr;
        bool             succ    = false;

        alloc.construct(newelem, e, currtop);

        PinGuard         guard(alloc, 1);

        // \mot
        // *1 Consume pin makes all stack updates made through the top
        //    visible to this thread.
        // *2 A successful c&e is a release operation (publish data)
        // *3 On failure, we load relaxed and fall-back to the consume load
        //
        // A model checked alternative would use
        //   (1) relaxed, (2) acq_rel, (3) relaxed
        // We posit that a consume fence would be a NOOP on most architectures,
        // and therefore cheaper.
        while (!succ)
        {
          stack_elem*    guardedptr = alloc.template pin<std::memory_order_consume>(top); // 1

          currtop = guardedptr;
          newelem->set_next(currtop);
          succ = top.compare_exchange_strong( currtop,
                                              newelem,
                                              std::memory_order_release,                  // 2
                                              std::memory_order_relaxed                   // 3
                                            );
          alloc.unpin(guardedptr, 0);
        }
      }

      std::pair<_Tp, bool> pop()
      {
        stack_elem*      nexttop = nullptr;
        stack_elem*      currtop = nullptr;
        bool             succ    = false;
        _Node_alloc_type alloc   = get_allocator();
        PinGuard         guard(alloc, 1);

        // \mot
        // *1 Consume pin makes all stack updates made through the top
        //    visible to this thread. The read of the currtop's next pointer
        //    relies on the consume fence.
        // *2 A successful c&e can be relaxed and we rely on a release sequence
        //    between the last push and any subsequent push and pop
        // *3 On failure, we load relaxed and fall-back to the consume pin
        //
        // A model checked alternative would use
        //   (1) relaxed, (2) acquire, (3) relaxed
        // We posit that a consume fence would be a NOOP on most architectures.
        while (!succ)
        {
          stack_elem*    guardedptr = alloc.template pin<std::memory_order_consume>(top); // 1

          if (guardedptr == nullptr) return std::make_pair(_Tp(), false);

          currtop = guardedptr;
          nexttop = currtop->next();
          succ = top.compare_exchange_strong( currtop,
                                              nexttop,
                                              std::memory_order_relaxed,                  // 2
                                              std::memory_order_relaxed                   // 3
                                            );
          alloc.unpin(guardedptr, 0);
        }

        ptr_guard<_Node_alloc_type> pguard(alloc, currtop);
        return std::make_pair(std::move(currtop->elem()), true);
      }

      _Node_alloc_type get_allocator() const { return nodeAlloc; }

    private:
      std::atomic<stack_elem*> top;
      _Node_alloc_type         nodeAlloc;
  };
}

namespace locking
{
  template <class T>
  struct stack
  {
      typedef aux::stack_elem<T>     stack_elem;
      typedef std::mutex             mutex;
      typedef std::lock_guard<mutex> lock_guard;

      stack()
      : m(), top(nullptr)
      {}

      void push(const T& e)
      {
        lock_guard  guard(m);
        stack_elem* currtop = top;

        top = new stack_elem(e, currtop);
      }

      std::pair<T, bool> pop()
      {
        lock_guard                  guard(m);
        std::unique_ptr<stack_elem> currtop(top);

        if (currtop == nullptr) return std::make_pair(T(), false);

        top = currtop->next();
        return std::make_pair(currtop->elem(), true);
      }

    private:
      mutex       m;
      stack_elem* top;
  };
}
