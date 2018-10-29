
/// \file    tasks.hpp
/// \brief   Framework for featherweight tasks.
/// \details Featherweight tasks are a portable yet effective prototype
///          framework for task scheduling and task reductions on shared memory
///          architectures. Featherweight tasks do not need their own stack
///          since they are executed on a regular thread's stack. Tasks may
///          spawn other tasks, or return a value contributing to a
///          reduction operation across all tasks. Currently, tasks cannot
///          wait for other tasks to be completed.
///
///          The framework has been designed with lock-free techniques and
///          generic programming principles in mind. Clients may plug-in a
///          variety of memory management techniques to customize the
///          framework for the application context ( see pmemory.hpp ).
///
///          Although, the featherweight tasks cannot wait on the completion of
///          other tasks, they can be applied in a variety of domains:
///          - reductions across tasks
///          - search space exploration
///          - partitioning of an iteration space
///
///          uab::execute_tasks() offers an example that implements task-based
///          computation of Fibonacci numbers.
///          Other examples can be found under BLAZE_HOME/examples/tasks .
/// \author  Peter Pirkelbauer ( pirkelbauer@uab.edu )

/*
 * A simple, portable, and generic framework for reductions over tasks.
 *
 * Implementer: Peter Pirkelbauer (UAB) - 2018
 *
 * This program is part of the Blaze Concurrent Library.
 * Copyright (c) 2018, University of Alabama at Birmingham
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation and/or
 * other materials provided with the distribution.
 * 3. Neither the name of the copyright holders nor the names of its contributors
 * may be used to endorse or promote products derived from this software without
 * specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 * OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 */


#ifndef _TASKS_HPP
#define _TASKS_HPP

#include <algorithm>
#include <thread>
#include <cassert>

// for setting and getting binding
#include <sched.h>
#include <pthread.h>

#include "atomicutil.hpp"
#include "pmemory.hpp"

#ifndef _ARCHMODEL_H
#include "archmodel.hpp"

typedef uab::generic_arch         arch_model;
#endif

namespace uab
{
  /// \brief void reduction defaults to empty operators
  /// \details
  ///    use:
  ///    \code
  ///    template <class Pool>
  ///    Void reductionOp(Pool taskpool, Task t)
  ///    {
  ///      return Void();
  ///    }
  ///    \endcode
  struct Void
  {
    Void() {}

    template <class T = int>
    Void(const T&) {}

    template <class T>
    Void operator+=(const T&) { return *this; }
  };

  //
  // Task pool implementation

  /// constant for block-size used in task pool
  static const size_t BLKSZ                = 1024;

  /// \private
  /// \brief Data block in a FIFO queue implemented as list of blocks.
  ///        dataQ supports a single producer, multiple consumers
  ///        tail->block<-block<-block<-head
  ///
  /// \details A data block holds BLKSZ elements at a time.
  ///          The data block is NOT cyclic, hence we get away with
  ///          relaxed memory operations. It is the task of the memory
  ///          management to guarantee that the deq operation completes
  ///          before the block can be recycled.
  /// \tparam T the task being stored
  template <class T>
  struct dataQ
  {
    /// type for head and tail position within block
    typedef uab::aligned_atomic_type<size_t, CACHELINESZ>    counter_t;

    /// pointer to next data block
    typedef uab::aligned_atomic_type<dataQ<T>*, CACHELINESZ> next_t;

    /// task type
    typedef T                                                value_t;

    static const size_t BLK = BLKSZ; ///< number of elements in a block

    counter_t  hd;         ///< block internal head
    counter_t  tl;         ///< block internal tail
    next_t     next;       ///< next block

    /// elements in the block
    /// \note in order to avoid interference from ctor and dtor calls
    ///       by classes with non-trivial constructor and deconstructor, the
    ///       block manages construction and deconstruction of its task objects.
    ///       Thus, it just reserves storage for BLK elements, and calls ctor
    ///       and dtor upon enqueue and dequeue (move based when possible).
    char       data[BLK*sizeof(value_t)];

    /// constructs data block with 0 elements
    dataQ()
    : hd(0), tl(0), next(nullptr)
    {}

    /// constructs a data block with 1 element
    /// \param el the first task
    explicit
    dataQ(value_t&& el)
    : hd(0), tl(1), next(nullptr)
    {
      new (at(0)) value_t (el);
    }

    /// constructs a data block with 1 element
    /// \param el the first task
    explicit
    dataQ(const value_t& el)
    : hd(0), tl(1), next(nullptr)
    {
      new (at(0)) value_t (el);
    }

    /// computes the address of element idx
    value_t* at(size_t idx)
    {
      char* res = data + idx*sizeof(T);

      assert(res >= data && res <= data + (BLKSZ-1)*sizeof(T));
      return reinterpret_cast<value_t*>(res);
    }

    /// tests if the data block can store more tasks
    bool has_space()
    {
      // \mo relaxed, b/c tl is only modified by this thread
      return tl.val.load(std::memory_order_relaxed) < BLK;
    }

    /// copy-enqueues a new element in this block
    /// \pre hasSpace() == true
    /// \pre only executed by the queue owner
    void enq(const value_t& el)
    {
      // \mo relaxed, b/c tl is only modified by this thread
      const size_t tail = tl.val.load(std::memory_order_relaxed);

      assert(tail < BLK);
      new (at(tail)) value_t (el); // in-place constructions

      // \mo release b/c data[tail] is published
      tl.val.store(tail+1, std::memory_order_release);
    }

    /// move-enqueues a new element in this block
    /// \pre hasSpace() == true
    /// \pre only executed by the queue owner
    void enq(value_t&& el)
    {
      // \mo relaxed, b/c tl is only modified by this thread
      const size_t tail = tl.val.load(std::memory_order_relaxed);

      assert(tail < BLK);
      new (at(tail)) value_t (el); // in-place constructions

      // \mo release b/c data[tail] is published
      tl.val.store(tail+1, std::memory_order_release);
    }

    /// dequeues a task and moves it into res
    /// \result the block from where the element was dequeued;
    ///         nullptr if unsuccessful.
    dataQ<T>* deq(value_t& res)
    {
      // \mo relaxed, b/c hd does not publish anything
      size_t head = hd.val.load(std::memory_order_relaxed);
      // \mo relaxed, b/c we are on the same thread
      size_t tail = tl.val.load(std::memory_order_relaxed);

      for (;;)
      {
        if (head == BLK)
        {
          // \mo consume, b/c we want to see initialized memory
          dataQ<T>* nxt = next.val.load(std::memory_order_relaxed);

          if (nxt == nullptr) return nullptr;

          return nxt->deq(res);
        }

        assert(head < BLK);
        if (head >= tail) return nullptr;

        // \mo relaxed (succ) since the data is never rewritten
        // \mo relaxed (fail) b/c nothing happened
        if (hd.val.compare_exchange_strong(head, head+1, std::memory_order_relaxed, std::memory_order_relaxed))
        {
          assert(head+1 <= tail);
          value_t* obj = at(head);

          // mv and deconstruct obj
          res = std::move(*obj);
          obj->~value_t();
          return this;
        }
      }
    }

    /// steals  a task and moves it into res
    /// \param  res the place where the task will be stored
    /// \param  tries number of attempts. Used to manage lock-free contention
    ///         by giving up after tries unsuccessful attempts to steal a task.
    /// \param  alloc the memory manager in use; needed to secure the next block
    ///         if traversal is needed.
    /// \result the block from where the element was dequeued;
    ///         nullptr if unsuccessful.
    /// \pre    this is pinned, if a finegrain (i.e., publishAndScan) is used
    /// \post   this is unpinned
    /// \note   publishAndScan is supported, though finegrain memory management
    ///         does not bring anything to the table for this task implementation.
    template <class Alloc>
    dataQ<T>* deq_steal(value_t& res, uint_fast8_t tries, Alloc alloc)
    {
      uint_fast8_t attempts = tries;

      // \mo relaxed, b/c hd does not publish anything
      size_t       head = hd.val.load(std::memory_order_relaxed);
      // \mo acquire, b/c we read from data[i], where i < tl
      size_t       tail = tl.val.load(std::memory_order_acquire);

      // give up under high contention
      while (attempts > 0)
      {
        if (head == BLK)
        {
          // \mo consume, b/c we want to see initialized memory
          dataQ<T>* nxt = alloc.template pin<std::memory_order_consume>(next.val);

          alloc.unpin(this, -1);

          if (nxt == nullptr) return nullptr;

          return nxt->deq_steal(res, tries, alloc);
        }

        assert(head < BLK);
        if (head >= tail) return nullptr;

        // \mo relaxed (succ) since the data is never rewritten
        // \mo relaxed (fail) b/c nothing happened
        if (hd.val.compare_exchange_strong(head, head+1, std::memory_order_relaxed, std::memory_order_relaxed))
        {
          assert(head+1 <= tail);
          value_t* obj = at(head);

          res = std::move(*obj);
          obj->~value_t();
          return this;
        }

        --attempts;
      }

      return nullptr;
    }
  };


  /// \private
  /// Main FIFO queue, built upon dataQ blocks + status data for threads
  /// \tparam T task type
  /// \tparam Alloc the memory manager as defined in pmemory.hpp
  /// \details
  ///    Each thread has its own flatQ object. The flatQ object also holds
  ///    some status information for its owning thread.
  template <class T, class Alloc>
  struct flatQ
  {
    /// Memory manager type
    typedef Alloc                                                        allocator_type;
    typedef std::allocator_traits<Alloc>                                 orig_alloc_traits;
    typedef typename orig_alloc_traits::template rebind_alloc<dataQ<T> > node_alloc_type;

    /// Number of concurrently pinned blocks (this and next)
    static const size_t                              SZPINWALL = 2;

    /// The allocator/memory manager
    node_alloc_type                                  nodealloc;

    /// keeps track of current
    bool                                             active;

    /// number of tasks stolen from this thread
    size_t                                           stolen;

    /// number of tasks enqueued by this thread
    size_t                                           tasks_made;

    /// number of completed tasks (incl. tasks solved by other threads, but
    ///   not tasks that this thread stole).
    size_t                                           tasks_done;

    /// tail (enqueue-side) of the queue; only accessed by owner
    dataQ<T>*                                        tail; // only accessed from queue owner

    /// head (dequeue-side) of the queue; concurrently accessed by work stealers.
    uab::aligned_atomic_type<dataQ<T>*, CACHELINESZ> head;

    /// estimate of available tasks
    uab::aligned_atomic_type<size_t, CACHELINESZ>    numtasks;

    /// number of stolen and completed tasks; updated by the work stealers.
    uab::aligned_atomic_type<size_t, CACHELINESZ>    steals;

    /// default ctor initializes an empty Q
    explicit
    flatQ(allocator_type alloc = allocator_type())
    : nodealloc(alloc), active(true), stolen(0), tasks_made(0), tasks_done(0),
      tail(new dataQ<T>()), head(tail), numtasks(0), steals(0)
    {}

    /// retrieve allocator
    node_alloc_type get_allocator()
    {
      return nodealloc;
    }

    /// allocates a new dataQ block
    dataQ<T>* new_block(T el)
    {
      node_alloc_type alloc = get_allocator();
      dataQ<T>*       tmp = alloc.allocate(1);

      alloc.construct(tmp, el);
      return tmp;
    }

    /// copy-enqueues a new element
    void enq(const T& el)
    {
      ++tasks_made;

      if (tail->has_space())
      {
        tail->enq(el);
        return;
      }

      dataQ<T>* tmp = new_block(el);

      numtasks.val.store(tasks_made - tasks_done, std::memory_order_relaxed);

      // \mo release, b/c new_block was initialized
      tail->next.val.store(tmp, std::memory_order_release);
      tail = tmp;
    }

    /// move-enqueues a new element
    void enq(T&& el)
    {
      ++tasks_made;

      if (tail->has_space())
      {
        tail->enq(el);
        return;
      }

      dataQ<T>* tmp = new_block(el);

      numtasks.val.store(tasks_made - tasks_done, std::memory_order_relaxed);

      // \mo release, b/c new_block was initialized
      tail->next.val.store(tmp, std::memory_order_release);
      tail = tmp;
    }

    /// records finished task
    void task_completed()
    {
      ++tasks_done;
    }

    /// tests if task became inactive (all tasks it produced (potentially stolen)
    ///   have been completed.
    bool became_inactive()
    {
      if (!active) return false;

      size_t prev_stolen = stolen;

      stolen = steals.val.load(std::memory_order_relaxed);

      assert(tasks_made - tasks_done >= stolen - prev_stolen);
      tasks_done += (stolen - prev_stolen);

      if (tasks_made > tasks_done) return false;

      active = false;
      return true;
    }

    /// dequeues an element and stores result into res
    /// \param  res  where the result will be stored
    /// \result true if successful, false otherwise
    bool deq(T& res)
    {
      typedef typename node_alloc_type::pinguard PinGuard;

      // \mo relaxed, b/c head is only updated by the owner
      dataQ<T>*               curr = head.val.load(std::memory_order_relaxed);
      assert(curr);

      // attempt to dequeue from own queue
      dataQ<T>* actl = curr->deq(res);

      // unsuccessful?
      if (actl == nullptr) return false;

      // all blocks that have become empty can be freed
      if (curr != actl)
      {
        node_alloc_type alloc = get_allocator();
        PinGuard        pinguard(alloc, SZPINWALL, lockfree::unordered());

        do
        {
          dataQ<T>* tmp = curr;

          curr = curr->next.val.load(std::memory_order_relaxed);
          alloc.deallocate(tmp, 1);
        } while (curr != actl);

        // \mo release to make content of new head available (needed?)
        head.val.store(actl, std::memory_order_release);

        // update available task counter
        numtasks.val.store(tasks_made - tasks_done, std::memory_order_relaxed);
      }

      return true;
    }

    /// attempts to steal a task from some other thread's queue
    /// \param  res   where the result will be stored
    /// \param  tries number of attempts before yielding to contention
    /// \result true if successful, false otherwise
    bool deq_steal(T& res, uint_fast8_t tries)
    {
      typedef typename node_alloc_type::pinguard PinGuard;

      node_alloc_type alloc = get_allocator();
      PinGuard        pinguard(alloc, SZPINWALL);

      // \mo consume, b/c we need to see memory initialized
      dataQ<T>*       curr = alloc.template pin<std::memory_order_consume>(head.val);
      assert(curr);

      return curr->deq_steal(res, tries, alloc) != nullptr;
    }

    /// quiescent release of memory
    /// \details
    ///    DO NOT CALL when concurrent operations are ongoing!
    void qrelease_memory()
    {
      get_allocator().qrelease_memory();
    }
  };


  /// The task pool class
  /// \tparam T task type
  /// \tparam _Alloc the memory manager (see pmemory.hpp)
  template <class T, class _Alloc = lockfree::epoch_manager<T, std::allocator>>
  struct pool
  {
    typedef T                                                       task_type;
    typedef _Alloc                                                  allocator_type;
    typedef flatQ<T, _Alloc>                                        tasq;

    typedef std::allocator_traits<_Alloc>                           orig_alloc_traits;
    typedef typename orig_alloc_traits::template rebind_alloc<tasq> node_alloc_type;

    /// maximum number of threads
    const uint_fast32_t                                  MAXTQ;

    /// the allocator/memory manager
    node_alloc_type                                      nodealloc;

    /// counts active threads; when count reaches 0 all task have been handled
    uab::aligned_atomic_type<uint_fast32_t, CACHELINESZ> active;

    /// task specific queue handles
    uab::aligned_atomic_type<tasq*, CACHELINESZ>         taskq[256]; // \todo remove magic constant

    /// thread id within pool
    static thread_local uint_fast32_t                    idx;

    /// victim (thread id) where the last task has been stolen
    static thread_local uint_fast32_t                    last_victim;

    /// this task queue
    static thread_local tasq*                            tq_loc;

    /// last victim's task queue
    static thread_local tasq*                            tq_rem;

    /// initializes pool and enqueues first task
    /// \param numthreads number of worker threads (and task queues)
    /// \param work       first task
    /// \param alloc      memory manager
    explicit
    pool(uint_fast32_t numthreads, T work, const allocator_type& alloc = allocator_type())
    : MAXTQ(numthreads), nodealloc(alloc), active(numthreads)
    {
      assert(tq_loc == nullptr);
      assert(MAXTQ <= 256);

      // initialize main thread
      tq_loc      = new tasq(nodealloc);
      idx         = 0;
      last_victim = 0;
      taskq[0].val.store(tq_loc, std::memory_order_relaxed);

      enq(std::move(work));

      // set tasqs to null
      for (uint_fast32_t i = 1; i < MAXTQ; ++i)
      {
        // \mo relaxed, b/c the pool is constructed before we fork
        taskq[i].val.store(nullptr, std::memory_order_relaxed);
      }
    }

    /// initializes thread local storage of task pool
    /// \param parent pool id of parent thread
    /// \param self   pool id of this thread
    /// \note pool ids are generated by the task system
    void init(uint_fast32_t parent, uint_fast32_t self)
    {
      if (self == 0) return;

      tq_loc      = new tasq(nodealloc);
      idx         = self;
      last_victim = parent;

      taskq[self].val.store(tq_loc, std::memory_order_release);
    }

    /// copy-enqueues a new task
    void enq(const T& el)
    {
      tq_loc->enq(el);
    }

    /// move-enqueues a new task
    void enq(T&& el)
    {
      tq_loc->enq(el);
    }

    /// called before a new task is scheduled
    void work_started()
    {
      if (!tq_loc->active)
      {
        active.val.fetch_add(1, std::memory_order_relaxed);
        tq_loc->active = true;
      }
    }

    /// called after a task has been completed
    void work_completed()
    {
      // was the task stolen?
      if (tq_rem == nullptr)
      {
        // no, just increment local work counter
        tq_loc->task_completed();
        return;
      }

      // yes, notify victim thread that its task has completed
      tq_rem->steals.val.fetch_add(1, std::memory_order_relaxed);
      tq_rem = nullptr;
    }

    /// checks whether all created tasks have completed
    /// \details
    ///    updates the pool's active counter if needed
    void check_active()
    {
      if (tq_loc->became_inactive())
      {
        active.val.fetch_sub(1, std::memory_order_relaxed);
      }
    }

    /// checks if there is at least one active thread remaining
    bool has_work()
    {
      return active.val.load(std::memory_order_relaxed) > 0;
    }

    /// dequeues a task and moves it into res
    /// \result  true if a task was dequeued, false otherwise
    bool deq(T& res)
    {
      // used to estimate number of stealing attempts
      static const size_t SAMPLESIZE = 4;

      // try to dequeue from the local queue first
      bool succ = tq_loc->deq(res);

      if (succ) return succ;

      // no more tasks in the local queue
      //   -> steal a task from somewhere else

      // estimate average number of tasks from a sample
      //   and choose first victim from the sample with the most tasks.
#if 1 /* informed stealing */
      size_t             tasks_sum = 0;
      size_t             tasks_avail[SAMPLESIZE];

      {
        uint_fast32_t    tasks_i = SAMPLESIZE;
        uint_fast32_t    probe[SAMPLESIZE] = { last_victim, last_victim+(idx+1), idx+(MAXTQ-1), idx+1 };
        size_t           max     = 0;

        while (tasks_i)
        {
          --tasks_i;
          uint_fast32_t  tid = probe[tasks_i];
          while (tid >= MAXTQ) { tid -= MAXTQ; }

          tasq*          victim = taskq[tid].val.load(std::memory_order_consume);

          if (victim)
          {
            tasks_avail[tasks_i] = victim->numtasks.val.load(std::memory_order_relaxed);
            tasks_sum += tasks_avail[tasks_i];

            if (tasks_avail[tasks_i] > max)
            {
              max = tasks_avail[tasks_i];
              last_victim = tid;
            }
          }
        }

        std::swap(tasks_avail[tasks_i], tasks_avail[0]);
      }
#endif /* informed stealing */
      {
        // work stealing
        uint_fast32_t      tasks_i = 0;
        uint_fast32_t      thrid   = last_victim;

        do
        {
          // \mo consume, b/c we follow the pointer
          tasq* victim = taskq[thrid].val.load(std::memory_order_consume);

          if (victim)
          {
            uint_fast8_t tries = arch_model::num_tries(idx, thrid);
#if 1 /* informed stealing */
            size_t        avail = victim->numtasks.val.load(std::memory_order_relaxed);

            if (avail * (SAMPLESIZE/2) > tasks_sum)
            {
              tries = std::min<uint_fast8_t>(tries+4, 8);
            }
            else if (avail * (SAMPLESIZE*4) < tasks_sum)
            {
              tries = 0; // 1, -2, /2, ...
            }

            tasks_sum -= tasks_avail[tasks_i];
            tasks_sum += (tasks_avail[tasks_i] = avail);
            ++tasks_i;
            if (tasks_i == SAMPLESIZE) tasks_i = 0;
#endif /* informed stealing */

            succ = victim->deq_steal(res, tries);

            if (succ)
            {
              tq_rem      = victim;
              last_victim = thrid;
              return succ;
            }
          }

          ++thrid;
          if (thrid == MAXTQ) thrid = 0;
        } while (thrid != last_victim);
      }

      last_victim = 0;
      assert(!succ);
      return succ;
    }

    /// releases all thread local memory
    void qrelease_memory()
    {
      tq_loc->qrelease_memory();
    }
  };

  template <class W, class A>
  thread_local
  uint_fast32_t pool<W,A>::idx = 0;

  template <class W, class A>
  thread_local
  uint_fast32_t pool<W,A>::last_victim = 0;

  template <class W, class A>
  thread_local
  typename pool<W,A>::tasq* pool<W,A>::tq_loc = nullptr;

  template <class W, class A>
  thread_local
  typename pool<W,A>::tasq* pool<W,A>::tq_rem = nullptr;

#if PRINT_STATS
  struct alignas(CACHELINESZ) threadstat
  {
    size_t num;        // number of tasks
    size_t core_init;  // core where thread began
    size_t core_last;  // core where thread ended
  };

  static threadstat work[NUMTHREADS];
#endif /* PRINT_STATS */

  /// \private
  /// \brief  main task loop
  /// \tparam R result type
  /// \tparam P pool type
  /// \tparam F function type
  template <class R, class P, class F>
  static
  void taskloop(size_t parent, size_t thrnum, R* res, P* tasks, F fun)
  {
    typedef typename P::task_type task_type;

    arch_model::bind_to_core(pthread_self(), thrnum);

#if PRINT_STATS
    work[thrnum].core_init = sched_getcpu();
#endif /* PRINT_STATS */

    assert(thrnum < tasks->MAXTQ);

    // initialize task-queues for each thread
    tasks->init(parent, thrnum);

    // initialize local reduction variable
    R      val = R();

    // run until no more tasks can be found
    do
    {
      task_type work;
      bool      succ = tasks->deq(work);

      while (succ)
      {
        tasks->work_started();
        val  += fun(*tasks, work);
        tasks->work_completed();

        succ  = tasks->deq(work);
      }

      tasks->check_active();
    } while (tasks->has_work());

    *res = val;

#if PRINT_STATS
    work[thrnum].num       = P::tq_loc->tasks_done;
    work[thrnum].core_last = sched_getcpu();
#endif

    // \todo consider detaching threads here
    tasks->qrelease_memory();
  }

  /// \private
  /// \brief  spawns n threads sequentially
  /// \tparam R result type
  /// \tparam P pool type
  /// \tparam F function type
  template <class R, class P, class F>
  static
  void spawn_threads_single(size_t parent, size_t thrnumlo, size_t thrnumhi, R* res, P* taskpool, F worker)
  {
    assert(thrnumhi > 0);
    const size_t last = thrnumhi-1;

    if (thrnumlo == last)
    {
      taskloop(parent, thrnumlo, res, taskpool, worker);
      return;
    }

    R sub;
    std::thread t(taskloop<R,P,F>, thrnumlo, last, &sub, taskpool, worker);

    spawn_threads_single<R>(thrnumlo, thrnumlo, last, res, taskpool, worker);
    t.join();

    *res += sub;
  }

  /// \private
  /// \brief  spawns threads; if the number exceeds a preset threshold
  ///         the thread spawns another spawning thread before spawning
  ///         half the threads.
  /// \tparam R result type
  /// \tparam P pool type
  /// \tparam F function type
  template <class R, class P, class F>
  static
  void spawn_threads(size_t parent, size_t thrnumlo, size_t thrnumhi, R* res, P* taskpool, F worker)
  {
    static const size_t SPAWN_THREASHOLD = 8;

    const size_t numthreads = thrnumhi - thrnumlo;

    if (numthreads <= SPAWN_THREASHOLD)
    {
      spawn_threads_single(parent, thrnumlo, thrnumhi, res, taskpool, worker);
      return;
    }

    const size_t thrnummid = thrnumlo+numthreads/2;
    R            sub;
    std::thread  spawner(spawn_threads<R,P,F>, thrnumlo, thrnummid, thrnumhi, &sub, taskpool, worker);

    spawn_threads<R>(parent, thrnumlo, thrnummid, res, taskpool, worker);
    spawner.join();

    *res += sub;
  }

/*
  /// \brief  spawns a spawning thread and continues to do work
  /// \tparam R result type
  /// \tparam P pool type
  /// \tparam F function type
  template <class R, class P, class F>
  static
  void spawn_spawner(size_t thrnumlo, size_t thrnumhi, R* res, P* taskpool, F worker)
  {
    static const size_t MIN_THREADS = 2;

    const size_t numthreads = thrnumhi - thrnumlo;

    if (numthreads < MIN_THREADS)
    {
      spawn_threads_single(thrnumlo, thrnumhi, res, taskpool, worker);
      return;
    }

    const size_t thrnummid = thrnumlo+1;
    R            sub;
    std::thread  spawner(spawn_threads<R,P,F>, thrnummid, thrnumhi, &sub, taskpool, worker);

    taskloop(thrnumlo, res, taskpool, worker);
    spawner.join();

    *res += sub;
  }
*/

  /// spawns numthread threads and calls fun(task)
  /// \tparam F type of functor
  /// \tparam T task type
  /// \param numthreads max number of threads to solve task
  /// \param fun functor that accepts a pool and a task
  /// \param task first task
  /// \details
  ///   example:
  /// \code
  ///    // fib functor
  ///    struct fib
  ///    {
  ///      // The task operator spawns fib(n-2) and continues with fib(n-1).
  ///      // The continuation-loop improves cache locality, and reduces
  ///      // overhead associated with enqueing and dequeing tasks.
  ///      template <class Pool>
  ///      int operator(Pool& p, int n)
  ///      {
  ///        if (n <= 1) return n;
  ///        while (--n > 1)
  ///        {
  ///          // spawn fib(n-2)
  ///          p.enq(n-1);
  ///
  ///          // continue with fib(n-1)
  ///        }
  ///        return 1; // reduced to base case
  ///      }
  ///    };
  ///
  ///    int res = execute_tasks(20 /* threads */, fib(), 10 /* argument to fib */);
  /// \endcode
  template <class F, class T>
  auto execute_tasks(size_t numthreads, F fun, T task) -> decltype(fun(*new pool<T>(0,task), task))
  {
    typedef decltype(fun(*new pool<T>(numthreads, task), task)) R;

    arch_model::set_threadinfo_info(numthreads);

    pool<T> taskpool(numthreads, task);
    R       res;

    spawn_threads(0, 0, numthreads, &res, &taskpool, fun);
    //~ spawn_spawner(0, numthreads, &res, &taskpool, fun);

#if PRINT_STATS
    for (size_t i = 0; i < NUMTHREADS; ++i)
    {
      std::cout << std::setw(2) << i << ": " << work[i].num
                << "   @" << work[i].core_init << " - " << work[i].core_last
                << std::endl;
    }
#endif

    return res;
  }

  /// spawns numthread threads and calls fun(task)
  /// \tparam A allocator/memory manager
  /// \tparam F type of functor
  /// \tparam T task type
  /// \param numthreads max number of threads to solve task
  /// \param fun functor that accepts a pool and a task
  /// \param task first task
  /// \details
  ///   example uses the gc_manager to handle deallocations by the
  ///   internal task pool.
  /// \code
  ///    int res = execute_tasks<gc_manager>(20, fib(), 10);
  /// \endcode
  template <class A, class F, class T>
  auto execute_tasks(size_t numthreads, F fun, T task) -> decltype(fun(*new pool<T,A>(0,task), task))
  {
    typedef decltype(fun(*new pool<T,A>(task), task)) R;

    arch_model::set_threadinfo_info(numthreads);

    pool<T> taskpool(task);
    R       res;

    spawn_threads(0, numthreads, &res, &taskpool, fun);

#if PRINT_STATS
    for (size_t i = 0; i < NUMTHREADS; ++i)
    {
      std::cout << std::setw(2) << i << ": " << work[i].num
                << "   @" << work[i].core_init << " - " << work[i].core_last
                << std::endl;
    }
#endif

    return res;
  }
} // end namespace uab

#endif /* _TASKS_HPP */
