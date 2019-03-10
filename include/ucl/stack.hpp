/// \file  stack.hpp
/// \brief Lock-free and lock-based stacks using linked lists as internal
///        representation.
/// \details
///        The implementation(s) rely on C++11.
///        - The lock-free stacks can be customized with memory managers ( see pmemory.hpp ).
///        - The lock-based stack can be customized with lock-type and guard types,
///          for example to elide the lock in the presence of hardware transactional
///          memory (requires HTM_ENABLED be set to 1).
/// \author Peter Pirkelbauer


#ifndef _UNLEASHED_STACK_HPP
#define _UNLEASHED_STACK_HPP 1

#include <atomic>
#include <mutex>
#include <utility>
#include <memory>

#include "pmemory.hpp"

/// \private
namespace aux
{
  /// \private
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
  /// \brief   Lock-free stack (Treiber) implementation for the relaxed memory model
  /// \details The stack is parametrized with available memory managers
  /// \tparam  _Tp the stack's value type
  /// \tparam  _Alloc the used memory manager (see pmemory.hpp )
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

      /// \brief constructs an empty stack with a given memory manager
      explicit
      stack(const allocator_type& alloc = allocator_type())
      : top(nullptr), nodeAlloc(alloc)
      {}

      /// \brief pushes a new element onto the stack
      void push(const _Tp& e)
      {
        _Node_alloc_type alloc   = get_allocator();
        stack_elem*      newelem = alloc.allocate(1);
        stack_elem*      currtop = nullptr;
        bool             succ    = false;

        alloc.construct(newelem, e, currtop);

        PinGuard         guard(alloc, 1);

        // \mo
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

      /// \brief pops an element from the stack
      /// \return a pair<_Tp, bool> where the flag indicates whether pop
      ///         was successful. If unsuccessful, first is default initialized.
      std::pair<_Tp, bool> pop()
      {
        stack_elem*      nexttop = nullptr;
        stack_elem*      currtop = nullptr;
        bool             succ    = false;
        _Node_alloc_type alloc   = get_allocator();
        PinGuard         guard(alloc, 1);

        // \mo
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

        ptr_deleter<_Node_alloc_type> elemguard(alloc, currtop);
        return std::make_pair(std::move(currtop->elem()), true);
      }

      /// \brief returns a copy of the allocator in use
      _Node_alloc_type
      get_allocator() const { return nodeAlloc; }

    private:
      std::atomic<stack_elem*> top;
      _Node_alloc_type         nodeAlloc;
  };
}

namespace locking
{
  /// \brief  a simple, lock-protected, linked-list based stack
  /// \tparam _Tp the stack's value type
  /// \tparam _M the used mutex class,
  /// \tparam _G the guard used, common options include std::lock_guard,
  ///            ucl::lockable_guard, ucl::elidable_guard. The chosen guard
  ///            needs to be compatible with the mutex interface.
  /// \details the stack uses new and delete to allocate its nodes.
  /// \todo    replace new/delete with allocators
  template <class _T, class _M = std::mutex, template <class> class _G = std::lock_guard>
  struct stack
  {
      typedef aux::stack_elem<_T> stack_elem;
      typedef _T                  value_type;
      typedef _M                  mutex;
      typedef _G<mutex>           lock_guard;

      /// \brief creates an empty stack
      stack()
      : m(), top(nullptr)
      {}

      /// \brief pushes a new element on to the stack
      void push(const value_type& e)
      {
        stack_elem* el = new stack_elem(e, nullptr);

        {
          lock_guard guard(m);

          el->set_next(top);
          top = el;
        }
      }

      /// \brief pops an element from the stack
      /// \return a pair<_Tp, bool> where the flag indicates whether pop
      ///         was successful. If unsuccessful, first is default initialized.
      std::pair<value_type, bool> pop()
      {
        stack_elem* el = nullptr;

        {
          lock_guard guard(m);
          el = top;

          if (el == nullptr) return std::make_pair(value_type(), false);
          top = el->next();
        }

        std::unique_ptr<stack_elem> oldptr(el);

        return std::make_pair(std::move(el->elem()), true);
      }

    private:
      mutex       m;
      stack_elem* top;
  };
}

#endif /* _UNLEASHED_STACK_HPP */
