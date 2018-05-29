#include <atomic>
#include <mutex>
#include <utility>
#include <cassert>

#include "pmemory.hpp"

namespace lockfree
{
  template <class _Tp>
  struct queue_elem
  {
      queue_elem(const _Tp& e)
      : elem(e), next(nullptr)
      {}

      _Tp                      elem;
      std::atomic<queue_elem*> next;
  };

  template < class _Tp,
             class _Alloc = lockfree::pub_scan_manager<_Tp, std::allocator>
           >
  struct queue
  {
      typedef queue_elem<_Tp> qelem_type;

      typedef typename _Alloc::value_type   _Alloc_value_type;
      typedef std::allocator_traits<_Alloc> _OrigAlloc_traits;
      typedef typename _OrigAlloc_traits::template rebind_alloc<qelem_type> _Node_alloc_type;
      typedef typename _Node_alloc_type::pinguard                           PinGuard;
      typedef _Alloc                                                        allocator_type;

      explicit
      queue(const allocator_type& alloc = allocator_type())
      : nodeAlloc(alloc), head(nodeAlloc.allocate(1)), tail(head.load(std::memory_order_relaxed))
      {}

      /// \brief  returns the allocator of this skiplist
      _Node_alloc_type get_allocator() const noexcept
      {
        return nodeAlloc;
      }      

      void enq(const _Tp& e)
      {
        _Node_alloc_type  alloc    = get_allocator();
        PinGuard          guard(alloc, 2);

        qelem_type* const entry    = new qelem_type(e);
        qelem_type*       last     = tail.load(std::memory_order_acquire);
        qelem_type*       expected = nullptr;

        while (!last->next.compare_exchange_strong(expected, entry, std::memory_order_release, std::memory_order_relaxed))
        {
          // Move tail to the new last element (help the other thread)
          //   if we fail some other thread must have succeeded before.
          const bool succ = tail.compare_exchange_strong(last, expected, std::memory_order_relaxed, std::memory_order_relaxed);

          // if we successfully replaced tail, we have a new last element
          //   otherwise the new last element was loaded by c&e
          if (succ) last = expected;

          // reset expected
          expected = nullptr;
        }

        // if we fail, some other thread must have updated tail meanwhile
        tail.compare_exchange_strong(last, entry, std::memory_order_release, std::memory_order_relaxed);
      }

      std::pair<_Tp, bool> deq()
      {
        qelem_type* dummy = head.load(std::memory_order_acquire);
        qelem_type* entry = nullptr;

        do
        {
          if (dummy == tail.load(std::memory_order_relaxed)) return std::make_pair(_Tp(), false);

          entry = dummy->next.load(std::memory_order_acquire);
          assert(entry);
        } while (!head.compare_exchange_strong(dummy, entry, std::memory_order_relaxed, std::memory_order_relaxed));

        return std::make_pair(entry->elem, true);
      }

    private:
      _Node_alloc_type         nodeAlloc;
      std::atomic<qelem_type*> head;
      std::atomic<qelem_type*> tail;
  };
}

namespace locking
{
  template <class T>
  struct queue_elem
  {
      queue_elem(const T& e)
      : elem(e), next(nullptr)
      {}

      T           elem;
      queue_elem* next;
  };

  template <class T>
  struct queue
  {
      typedef std::mutex             mutex;
      typedef std::lock_guard<mutex> lock_guard;

      queue()
      : m(), head(new queue_elem<T>(T())), tail(head)
      {}

      void enq(const T& e)
      {
        lock_guard     guard(m);
        queue_elem<T>* const entry = new queue_elem<T>(e);

        tail->next = entry;
        tail = entry;
      }

      std::pair<T, bool> deq()
      {
        lock_guard  guard(m);

        if (head == tail) return std::make_pair(T(), false);

        std::unique_ptr<queue_elem<T> > ptrdel( head );

        head = head->next;
        return std::make_pair(head->elem, true);
      }

    private:
      mutex          m;
      queue_elem<T>* head;
      queue_elem<T>* tail;
  };
}
