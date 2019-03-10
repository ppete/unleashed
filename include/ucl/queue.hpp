/// \file  queue.hpp
/// \brief Lock-free and lock-based queues using linked lists as internal
///        representation.
/// \details
///        The implementation(s) rely on C++11.
///        - The lock-free queue can be customized with memory managers ( see pmemory.hpp ).
///        - The lock-based queue can be customized with lock-type and guard types,
///          for example to elide the lock in the presence of hardware transactional
///          memory (requires HTM_ENABLED be set to 1).
/// \author Peter Pirkelbauer

#ifndef _UNLEASHED_QUEUE_HPP
#define _UNLEASHED_QUEUE_HPP 1

#include <atomic>
#include <mutex>
#include <utility>
#include <cassert>

#include "pmemory.hpp"

namespace lockfree
{
  /// \private
  template <class _Tp>
  struct queue_elem
  {
      queue_elem(_Tp&& e)
      : elem(std::move(e)), next(nullptr)
      {}

      _Tp                      elem;
      std::atomic<queue_elem*> next;
  };


  /// \brief   Lock-free FIFO queue (Michael/Scott) implemented for the C++11
  ///          relaxed memory model.
  /// \details The queue is parametrized with available memory managers
  /// \tparam  _Tp the stack's value type
  /// \tparam  _Alloc the used memory manager (see pmemory.hpp )
  template < class _Tp,
             class _Alloc = lockfree::pub_scan_manager<_Tp, std::allocator>
           >
  struct queue
  {
      typedef _Tp             value_type;
      typedef queue_elem<_Tp> qelem_type;

      typedef typename _Alloc::value_type   _Alloc_value_type;
      typedef std::allocator_traits<_Alloc> _OrigAlloc_traits;
      typedef typename _OrigAlloc_traits::template rebind_alloc<qelem_type> _Node_alloc_type;
      typedef typename _Node_alloc_type::pinguard                           PinGuard;
      typedef _Alloc                                                        allocator_type;

      /// \brief constructs an empty queue with a given memory manager
      explicit
      queue(const allocator_type& alloc = allocator_type())
      : nodeAlloc(alloc), head(nodeAlloc.allocate(1)), tail(head.load(std::memory_order_relaxed))
      {
        new (head.load(std::memory_order_relaxed)) qelem_type (_Tp());
      }

      /// \brief  returns the allocator of this skiplist
      _Node_alloc_type get_allocator() const noexcept
      {
        return nodeAlloc;
      }

      /// \brief enqueues a new element
      void enq(value_type&& e)
      {
        _Node_alloc_type  alloc    = get_allocator();
        PinGuard          guard(alloc, 2);
        qelem_type* const newelem  = alloc.allocate(1);

        //~ alloc.construct(newelem, std::move(e));
        new (newelem) qelem_type (std::move(e));

        qelem_type*       last     = alloc.template pin<std::memory_order_relaxed>(tail);
        qelem_type*       expected = nullptr;

        while (!last->next.compare_exchange_strong( expected,
                                                    newelem,
                                                    std::memory_order_acq_rel,
                                                    std::memory_order_acq_rel
                                                   ))
        {
          // Move tail to the new last element (help the other thread)
          //   if we fail some other thread must have succeeded before.
          const bool succ = tail.compare_exchange_strong( last,
                                                          expected,
                                                          std::memory_order_release,
                                                          std::memory_order_relaxed
                                                        );

          if (is_alloc_kind<finegrain>(typename _Node_alloc_type::manager_kind()))
          {
            alloc.unpin_all();

            // \todo this should not be a complete pin, but just a confirm
            last = alloc.template pin<std::memory_order_relaxed>(tail);
          }
          else
          {
            // if we successfully replaced tail, we have a new last element
            //   otherwise the new last element was loaded by c&e
            if (succ) last = expected;
          }

          // reset expected
          expected = nullptr;
        }

        // if we fail, some other thread must have updated tail meanwhile
        tail.compare_exchange_strong( last,
                                      newelem,
                                      std::memory_order_release,
                                      std::memory_order_relaxed
                                    );
      }

      /// \brief dequeues an element
      /// \return a pair<_Tp, bool> where the flag indicates whether the dequeue
      ///         was successful. If unsuccessful, first is default initialized.
      std::pair<_Tp, bool> deq()
      {
        _Node_alloc_type alloc = get_allocator();
        PinGuard         guard(alloc, 2);
        qelem_type*      dummy = nullptr;
        qelem_type*      entry = nullptr;

        do
        {
          alloc.unpin_all();

          dummy = alloc.template pin<std::memory_order_acquire>(head);
          if (dummy == tail.load(std::memory_order_acquire))
          {
            return std::make_pair(_Tp(), false);
          }

          entry = alloc.template pin<std::memory_order_relaxed>(dummy->next);
          assert(entry);
        } while (!head.compare_exchange_strong( dummy,
                                                entry,
                                                std::memory_order_release,
                                                std::memory_order_relaxed
                                               ));

        alloc.deallocate(dummy, 1);
        return std::make_pair(std::move(entry->elem), true);
      }

      void qrelease_memory()
      {
        get_allocator().qrelease_memory();
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
      queue_elem(T&& e)
      : elem(std::move(e)), next(nullptr)
      {}

      T           elem;
      queue_elem* next;
  };

  /// \brief  a simple, single lock-protected, linked-list based queue
  /// \tparam _T the stack's value type
  /// \tparam _M the used mutex class,
  /// \tparam _G the guard used, common options include std::lock_guard,
  ///            ucl::lockable_guard, ucl::elidable_guard. The chosen guard
  ///            needs to be compatible with the mutex interface.
  /// \details the stack uses new and delete to allocate its nodes.
  /// \todo    replace new/delete with allocators
  template <class _T, class _M = std::mutex, template <class> class _G = std::lock_guard>
  struct queue
  {
      typedef _T                                      value_type;
      typedef _M                                      mutex;
      typedef _G<mutex>                               lock_guard;
      typedef std::allocator<queue_elem<value_type> > allocator_type;

      queue()
      : m(), head(new queue_elem<value_type>(value_type())), tail(head)
      {}

      queue(allocator_type)
      : queue()
      {}

      void enq(value_type&& e)
      {
        queue_elem<value_type>* const entry = new queue_elem<value_type>(std::move(e));

        {
          lock_guard     guard(m);

          tail->next = entry;
          tail = entry;
        }
      }

      std::pair<value_type, bool> deq()
      {
        lock_guard  guard(m);

        if (head == tail) return std::make_pair(value_type(), false);

        std::unique_ptr<queue_elem<value_type> > ptrdel( head );

        head = head->next;
        return std::make_pair(std::move(head->elem), true);
      }

      /// for compatibility
      void qrelease_memory() {}

    private:
      mutex                   m;
      queue_elem<value_type>* head;
      queue_elem<value_type>* tail;
  };
}

#endif /* _UNLEASHED_QUEUE_HPP */
