
/// \file    task.hpp
/// \brief   Framework for featherweight tasks.
/// \details Featherweight tasks are a portable yet effective prototype
///          framework for task scheduling and task reductions on shared memory
///          architectures. Featherweight tasks do not need their own stack
///          since they are executed on a regular thread's stack. Tasks may
///          spawn other tasks, or return a value contributing to a
///          reduction operation across all tasks.
///
///          The framework has been designed with lock-free techniques and
///          generic programming principles in mind.
///          - Clients may plugin different task pools.
///          - Task pools, in turn, can be customized with different memory
///            reclamation techniques in order to choose a suitable mechanism
///            for different environments ( see pmemory.hpp ).
///
///          Although, the featherweight tasks cannot wait on the completion of
///          other tasks, they can be applied in a variety of domains:
///          - reductions across tasks
///          - search space exploration
///          - partitioning of an iteration space
///
///          If waiting is necessary, users may use continuations
///          (e.g., continuation or xcontinuation) to achieve coordination
///          among tasks.
///
///          ucl::execute_tasks() offers an example that implements
///          task-based computation of Fibonacci numbers.
///          Other examples can be found under UCL_HOME/examples/tasks .
/// defines Concept definitions for the task library:
/// \code{.cpp}
///
/// // defines the TaskFunctor concept
/// concept TaskFunctor
/// {
///   // result type of a task's execution
///   // \note the type is inferred and a type with that name does not
///   //       need to be present.
///   typename result_type;
///
///   // the executer of a task, takes a pool and its value_type
///   template <class TaskPool>
///   result_type operator()(TaskPool& pool, TaskPool::value_type&& task);
/// }
///
/// // Single producer, multiple consumer data structure.
/// //   The data structure is unordered in principle.
/// concept SPMCBag
/// {
///   // defines the task type this pool is able to handle
///   typename value_type;
///
///   // defines the allocator type in use
///   typename allocator_type;
///
///   // constructor that takes an allocator
///   SPMCBag(allocator_type alloc);
///
///   // adds a new task t to the bag
///   void enq(value_type t);
///
///   // dequeues a task and moves it into res
///   // \result true if a task was dequeued, false otherwise
///   bool deq(T& res);
///
///   // returns an estimated number of available tasks
///   size_t estimate_size();
///
///   // does any created task (stolen or not) remain unfinished?
///   bool has_unfinished_tasks();
///
///   // records finished task
///   void task_completed();
///
///   // records finished stolen task
///   void stolen_task_completed();
///
///   // releases resources after all tasks have been handled
///   // \details
///   //    it can be assumed that this method is called during quiescent time
///   void qrelease_memory();
/// }
///
/// // The taskpool concept
/// concept TaskPool
/// {
///   // defines the task type this pool is able to handle
///   typename value_type;
///
///   // returns the number of worker threads
///   size_t num_threads();
///
///   // initializes pool thread local data for thisthread
///   // and sets its parent thread
///   void init(size_t parent, size_t thisthread);
///
///   // attempts to get a task from some queue and stores it in t
///   // \result true, iff a task could be retrieved
///   bool deq(value_type& t);
///
///   // adds a new task t to this thread's queue
///   void enq(value_type t);
///
///   // called before work on task starts
///   void work_started();
///
///   // called after work on task finished
///   void work_completed();
///
///   // returms true if there is at least one active task in the system
///   bool work_available();
///
///   // releases (thread local) resources after all tasks have been handled
///   // \details
///   //    it can be assumed that this method is called during quiescent time
///   void qrelease_memory();
/// }
/// \endcode
///
/// \author  Peter Pirkelbauer

/*
 * A simple, portable, and generic framework for reductions over tasks.
 *
 * Implementer: Peter Pirkelbauer - 2018
 *
 * This program is part of the Unleashed Concurrency Library (UCL).
 *
 * Copyright (c) 2019, Peter Pirkelbauer
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


#ifndef _UCL_TASK_HPP
#define _UCL_TASK_HPP 1

#include <algorithm>
#include <thread>
#include <cassert>
#include <functional>

// for setting and getting binding
#include <sched.h>
#include <pthread.h>

#include "atomicutil.hpp"
#include "pmemory.hpp"

#ifndef _ARCHMODEL_H
#include "archmodel.hpp"

typedef ucl::generic_arch         arch_model;
#endif

// Switches on (1) and off (0) code for estimating the number of
// elements in task pools.
#define INFORMED_STEALING 1



namespace ucl
{
  /// \brief void reduction defaults to empty operators
  /// \details
  ///    use:
  ///    \code{.cpp}
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
  static const size_t BLKSZ                = 512;

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
    typedef ucl::aligned_atomic_type<size_t>    counter_t;

    /// pointer to next data block
    typedef ucl::aligned_atomic_type<dataQ<T>*> next_t;

    /// task type
    typedef T                                   value_t;

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
      new (at(0)) value_t (std::move(el));
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
      new (at(tail)) value_t (std::move(el)); // in-place constructions

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
  template <class T, class UclAlloc>
  struct flatQ
  {
      /// Memory manager type
      typedef UclAlloc                                                     orig_allocator_type;
      typedef std::allocator_traits<UclAlloc>                              orig_alloc_traits;
      typedef typename orig_alloc_traits::template rebind_alloc<dataQ<T> > allocator_type;
      typedef T                                                            value_type;

      /// Number of concurrently pinned blocks (this and next)
      static const size_t                              SZPINWALL = 2;

    // private:
      /// The allocator/memory manager
      allocator_type                                   nodealloc;

      /// tail (enqueue-side) of the queue; only accessed by owner
      dataQ<T>*                                        tail; // only accessed from queue owner

      /// head (dequeue-side) of the queue; concurrently accessed by work stealers.
      ucl::aligned_atomic_type<dataQ<T>*, CACHELINESZ> head;

      /// number of tasks stolen from this thread
      size_t                                           stolen;

      /// number of tasks enqueued by this thread
      size_t                                           tasks_made;

      /// number of completed tasks (incl. tasks solved by other threads, but
      ///   not tasks that this thread stole).
      size_t                                           tasks_done;

      /// number of stolen and completed tasks; updated by the work stealers.
      ucl::aligned_atomic_type<size_t, CACHELINESZ>    steals;

      /// estimate of available tasks
      ucl::aligned_atomic_type<size_t, CACHELINESZ>    numtasks;

    public:
      /// default ctor initializes an empty Q
      explicit
      flatQ(orig_allocator_type alloc = orig_allocator_type())
      : nodealloc(alloc), tail(nullptr), head(nullptr), stolen(0),
        tasks_made(0), tasks_done(0), steals(0), numtasks(0)
      {
        nodealloc.initialize_if_needed();

        tail = nodealloc.allocate(1);
        new (tail) dataQ<T>();
        head.val.store(tail);
      }

      /// retrieve allocator
      allocator_type get_allocator()
      {
        return nodealloc;
      }

      /// returns an estimated number of available tasks
      size_t estimate_size() const
      {
        return numtasks.val.load(std::memory_order_relaxed);
      }

      /// allocates a new dataQ block
      dataQ<T>* new_block(T&& el)
      {
        allocator_type  alloc = get_allocator();
        dataQ<T>*       tmp = alloc.allocate(1);

        new (tmp)dataQ<T>(std::move(el));
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

        dataQ<T>* tmp = new_block(std::move(T(el)));

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
          tail->enq(std::move(el));
          return;
        }

        dataQ<T>* tmp = new_block(std::move(el));

        numtasks.val.store(tasks_made - tasks_done, std::memory_order_relaxed);

        // \mo release, b/c new_block was initialized
        tail->next.val.store(tmp, std::memory_order_release);
        tail = tmp;
      }

      /// dequeues an element and stores result into res
      /// \param  res  where the result will be stored
      /// \result true if successful, false otherwise
      bool deq(T& res)
      {
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
          allocator_type alloc = get_allocator();

          do
          {
            dataQ<T>* tmp = curr;

            curr = curr->next.val.load(std::memory_order_relaxed);
            alloc.deallocate(tmp, 1);
          } while (curr != actl);

          // \mo release to make content of new head available (needed?)
          head.val.store(actl, std::memory_order_release);

#if INFORMED_STEALING
          // update completed stolen tasks
          has_unfinished_tasks();
#endif /* INFORMED_STEALING */

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
        // \mo consume, b/c we need to see memory initialized
        allocator_type  alloc = get_allocator();
        dataQ<T>*       curr  = alloc.template pin<std::memory_order_consume>(head.val);
        assert(curr);

        return curr->deq_steal(res, tries, alloc) != nullptr;
      }

      /// records finished task
      void task_completed() { ++tasks_done; }

      /// records finished stolen task
      void stolen_task_completed()
      {
        steals.val.fetch_add(1, std::memory_order_relaxed);
      }

      /// does any created task (stolen or not) remain unfinished?
      bool has_unfinished_tasks();

      /// quiescent release of memory
      /// \details
      ///    DO NOT CALL when concurrent operations are ongoing!
      void qrelease_memory()
      {
        get_allocator().qrelease_memory();
      }
  };

  template <class T, class UclAlloc>
  bool flatQ<T, UclAlloc>::has_unfinished_tasks()
  {
    size_t prev_stolen = stolen;

    stolen = steals.val.load(std::memory_order_relaxed);

    assert(tasks_made - tasks_done >= stolen - prev_stolen);
    tasks_done += (stolen - prev_stolen);

    return tasks_made > tasks_done;
  }


  template <class Q>
  struct active_tracker : Q
  {
    /// keeps track of current
    bool                    active;

#if UCL_RUNTIME_DATA
    size_t                  task_acq = 0;
    size_t                  task_acq_tries = 0;
#endif /* UCL_RUNTIME_DATA */

    active_tracker()
    : Q(), active(true)
    {}

    template <class UclAlloc>
    active_tracker(UclAlloc alloc)
    : Q(alloc), active(true)
    {}
  };


  /// The task pool class
  /// \tparam SPMCBag a single-producer multi-consumer data structure
  ///         for managing tasks.
  ///         The owning thread is the single producer, while all
  ///         threads may be consumers.
  template <class SPMCBag>
  struct perthread_pool
  {
    typedef active_tracker<SPMCBag>                      tasq;
    typedef typename SPMCBag::value_type                 value_type;
    typedef typename SPMCBag::allocator_type             allocator_type;

    /// maximum number of threads
    const uint_fast32_t                                  MAXTQ;

    /// the allocator/memory manager
    allocator_type                                       nodealloc;

    /// counts active threads; when count reaches 0 all task have been handled
    ucl::aligned_atomic_type<uint_fast32_t, CACHELINESZ> active;

    /// task specific queue handles
    ucl::aligned_atomic_type<tasq*, CACHELINESZ> taskq[256]; // \todo remove magic constant

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
    perthread_pool(uint_fast32_t numthreads, value_type work, const allocator_type& alloc = allocator_type())
    : MAXTQ(numthreads), nodealloc(alloc), active(numthreads)
    {
      assert(tq_loc == nullptr);
      assert(MAXTQ <= 256);

      // initialize main thread
      // \todo \new
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

    /// retrieve allocator
    allocator_type get_allocator()
    {
      return nodealloc;
    }

    /// returns the maximum number of worker threads for this pool
    size_t num_threads() const { return MAXTQ; }

    /// initializes thread local storage of task pool
    /// \param parent pool id of parent thread
    /// \param self   pool id of this thread
    /// \note pool ids are generated by the task system
    void init(uint_fast32_t parent, uint_fast32_t self)
    {
      nodealloc.initialize_if_needed();

      if (self == 0) return;

      // \todo \new
      tq_loc      = new tasq(nodealloc);
      idx         = self;
      last_victim = parent;

      taskq[self].val.store(tq_loc, std::memory_order_release);
    }

    /// copy-enqueues a new task
    void enq(const value_type& el)
    {
      tq_loc->enq(el);
    }

    /// move-enqueues a new task
    void enq(value_type&& el)
    {
      tq_loc->enq(std::move(el));
    }

    /// called before a new task is scheduled
    void work_started()
    {
      if (tq_loc->active) return;

      active.val.fetch_add(1, std::memory_order_relaxed);
      tq_loc->active = true;
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
      tq_rem->stolen_task_completed();
      tq_rem = nullptr;
    }

    bool became_inactive(tasq* q)
    {
      if (!q->active || q->has_unfinished_tasks())
        return false;

      q->active = false;
      return true;
    }

    /// checks if there is at least one active thread remaining
    bool work_available()
    {
      size_t cntactv = became_inactive(tq_loc)
                          ? active.val.fetch_sub(1, std::memory_order_relaxed) - 1
                          : active.val.load(std::memory_order_relaxed)
                          ;

      return cntactv > 0;
    }

    /// dequeues a task and moves it into res
    /// \result  true if a task was dequeued, false otherwise
    bool deq(value_type& res)
    {
      // used to estimate number of stealing attempts

      // try to dequeue from the local queue first
      bool succ = tq_loc->deq(res);

      if (succ) return succ;

      // no more tasks in the local queue
      //   -> steal a task from somewhere else

#if INFORMED_STEALING
      // estimate average number of tasks from a sample
      //   and choose first victim from the sample with the most tasks.
      static const size_t SAMPLESIZE = 4;
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
            tasks_avail[tasks_i] = victim->estimate_size();
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
#endif /* INFORMED_STEALING */
      {
        typedef typename allocator_type::pinguard PinGuard;

        allocator_type     alloc = get_allocator();
        PinGuard           pinguard(alloc, SPMCBag::SZPINWALL);

        // work stealing
        uint_fast32_t      thrid   = last_victim;

#if INFORMED_STEALING
        uint_fast32_t      tasks_i = 0;
#endif /* INFORMED_STEALING */

        do
        {
          // \mo consume, b/c we follow the pointer
          tasq* victim = taskq[thrid].val.load(std::memory_order_consume);

          if (victim)
          {
            uint_fast8_t tries = arch_model::num_tries(idx, thrid);

#if INFORMED_STEALING
            size_t        avail = victim->estimate_size();

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
#endif /* INFORMED_STEALING */

            succ = victim->deq_steal(res, tries);

            if (succ)
            {
              assert(idx != thrid);
#if UCL_RUNTIME_DATA
              ++tq_loc->task_acq;
#endif /* UCL_RUNTIME_DATA */

              tq_rem      = victim;
              last_victim = thrid;
              return succ;
            }

#if UCL_RUNTIME_DATA
              tq_loc->task_acq_tries += tries;
#endif /* UCL_RUNTIME_DATA */
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

  template <class Q>
  thread_local
  uint_fast32_t perthread_pool<Q>::idx = 0;

  template <class Q>
  thread_local
  uint_fast32_t perthread_pool<Q>::last_victim = 0;

  template <class Q>
  thread_local
  typename perthread_pool<Q>::tasq* perthread_pool<Q>::tq_loc = nullptr;

  template <class Q>
  thread_local
  typename perthread_pool<Q>::tasq* perthread_pool<Q>::tq_rem = nullptr;


  /// \private
  /// \brief  main task loop
  /// \tparam R result type
  /// \tparam P pool type
  /// \tparam F function type
  template <class R, class P, class F>
  static
  void taskloop(size_t parent, size_t thrnum, R* res, P* tasks, F fun)
  {
    typedef typename P::value_type task_type;

#if !defined __CYGWIN__ && !defined OS_WIN32 && !defined __OpenBSD__ && !defined __sun
    arch_model::bind_to_core(pthread_self(), thrnum);
#endif

    assert(thrnum < tasks->num_threads());

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
        val  += fun(*tasks, std::move(work));
        tasks->work_completed();

        succ  = tasks->deq(work);
      }
    } while (tasks->work_available());

    *res = val;

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

#if UCL_RUNTIME_DATA
  template <class T>
  static inline
  void log_telemetry(std::ostream& os, ucl::perthread_pool<T>& pool)
  {
    typedef typename ucl::perthread_pool<T>::tasq tasq;

    size_t total_stolen = 0;
    size_t total_tasks  = 0;
    size_t total_task_acq = 0;
    size_t total_task_acq_tries = 0;

    for (size_t i = 0; i < pool.MAXTQ; ++i)
    {
      tasq& tsk = *pool.taskq[i].val;

      total_tasks          += tsk.tasks_done;
      total_stolen         += tsk.stolen;
      total_task_acq       += tsk.task_acq;
      total_task_acq_tries += tsk.task_acq_tries;

      os << "* t" << i
         << "   " << tsk.tasks_done
         << " (" << tsk.stolen << "/" << tsk.task_acq << "/" << tsk.task_acq_tries
         << ")"
         << std::endl;
    }

    os << "**** totals: " << total_tasks
       << " (" << total_stolen << " " << ((total_stolen * 100.0) / (total_tasks))
       << "%/" << total_task_acq << "/" << total_task_acq_tries << ")"
       << std::endl;
  }
#endif /* UCL_RUNTIME_DATA */


  /// executes fun over the taskpool
  /// \tparam TaskFunctor a functor being able to execute tasks in a TaskPool
  ///         (required to implement the TaskFunctor concept).
  /// \tparam TaskPool a class that keeps track of available tasks
  ///         (required to implement the TaskPool concept).
  /// \param fun functor that accepts a pool and a task
  /// \param pool a task-pool
  template <class TaskFunctor, class TaskPool>
  auto execute_pool(TaskFunctor fun, TaskPool& pool)
       -> decltype( fun(pool, std::move(typename TaskPool::value_type())) )
  {
    typedef decltype( fun(pool, std::move(typename TaskPool::value_type())) ) R;

    size_t numthreads = pool.num_threads();
    R      res;

    arch_model::set_threadinfo_info(numthreads);

    spawn_threads(0, 0, numthreads, &res, &pool, fun);
#if UCL_RUNTIME_DATA
    log_telemetry(std::cerr, pool);

    log_epoch_telemetry(std::cerr, pool.get_allocator());
#endif /* UCL_RUNTIME_DATA */
    return res;
  }

  template <class T>
  using default_task_allocator = lockfree::epoch_manager<T, std::allocator>;

  template <class T>
  using fifo_queue_pool = ucl::perthread_pool<flatQ<T, default_task_allocator<T> > >;

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
  template <class F, class T, class UclAlloc = lockfree::epoch_manager<T, std::allocator> >
  static inline
  auto execute_tasks(size_t numthreads, F fun, T&& task)
       -> decltype(execute_pool(fun, *new fifo_queue_pool<T>(numthreads, std::move(task))))
  {
    fifo_queue_pool<T> pool(numthreads, std::move(task));

    return execute_pool(fun, pool);
  }



#if OBSOLETE
  use execute_tasks above

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
  static inline
  auto execute_tasks(size_t numthreads, F fun, T&& task)
       -> decltype(execute_pool(fun, *new perthread_pool<T,A>(numthreads, std::move(task))))
  {
    perthread_pool<T,A> taskpool(numthreads, std::move(task));

    return execute_pool(fun, taskpool);
  }
#endif

} // end namespace ucl

#endif /* _UCL_TASK_HPP */
