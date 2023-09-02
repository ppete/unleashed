/*
 * A simple, portable, and generic framework for reductions over tasks.
 *
 * Implementer: Peter Pirkelbauer
 *
 * This program is part of the Unleashed Concurrency Library (UCL).
 *
 * Copyright (c) 2019-2020, Peter Pirkelbauer
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
/// defines concept signatures for the task library:
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
///   bool deq(void* rawmem);
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
///   // \param resetGlobals if set, the thread will be responsible for
///   //                     clearing all global state. Only one thread should receive
///   //                     resetGlobals.
///   // \details
///   //    it can be assumed that this method is called during quiescent time
///   void qshutdown(bool resetGlobals);
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
///   bool deq(void* rawmem);
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
///   // \param resetGlobals if set, the thread will be responsible for
///   //                     clearing all global state. Only one thread should receive
///   //                     resetGlobals.
///   // \details
///   //    it can be assumed that this method is called during quiescent time
///   void qshutdown(bool resetGlobals);
/// }
/// \endcode
///
/// \author  Peter Pirkelbauer

#ifndef _UCL_TASK_HPP
#define _UCL_TASK_HPP 1

#define _GNU_SOURCE 1

#include <algorithm>
#include <cassert>
#include <functional>

#if UCL_RUNTIME_DATA
#include <sys/time.h>
#include <sys/resource.h>
#endif /* UCL_RUNTIME_DATA */

#include "ucl/thread.hpp"
#include "ucl/atomicutil.hpp"
#include "ucl/pmemory.hpp"

#ifndef _UCL_ARCHMODEL_H
#include "ucl/archmodel.hpp"

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
  
  //~ typedef int nat_t;          /* int faster than size_t? */
  //~ typedef int_fast64_t nat_t; /* int faster than size_t? */
  //~ typedef int_fast32_t nat_t; /* int faster than size_t? */
  typedef size_t nat_t;

  /// constant for default block-size (number of elements) in the task pool.
  static constexpr nat_t DEFAULT_BLOCK_SIZE = 1024;

#if UCL_RUNTIME_DATA
  struct ucl_runtime_data
  {
    nat_t minor_page_faults = 0;
    
    ucl_runtime_data& 
    operator+=(const ucl_runtime_data& other)
    {
      this->minor_page_faults += other.minor_page_faults;
      return *this;
    }
  };
  
  std::ostream& operator<<(std::ostream& os, ucl_runtime_data data)
  {
    return os << data.minor_page_faults;
  }

  thread_local ucl_runtime_data tk_enq_data;
  thread_local ucl_runtime_data tk_deq_steal_data;
  thread_local ucl_runtime_data tk_deq_data;
  thread_local ucl_runtime_data tk_deq_collect_data;

  thread_local nat_t task_steals      = 0;
  thread_local nat_t task_steal_tries = 0;
  thread_local nat_t task_steal_skips = 0;
  
  void ucl_rt_begin(ucl_runtime_data& rec)
  {
    rusage usage;
    
    getrusage(RUSAGE_THREAD, &usage);
    
    rec.minor_page_faults = usage.ru_minflt;
  }
  
  void ucl_rt_end(ucl_runtime_data& collect, const ucl_runtime_data& rec)
  {
    rusage usage;
    
    getrusage(RUSAGE_THREAD, &usage);
    
    collect.minor_page_faults += usage.ru_minflt;
    collect.minor_page_faults -= rec.minor_page_faults;
  }
#endif /* UCL_RUNTIME_DATA */

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
    typedef ucl::aligned_atomic_type<nat_t>           counter_t;

    /// pointer to next data block
    typedef ucl::aligned_atomic_type<dataQ<T>*>       next_t;

    /// task type
    typedef T                                         value_t;

    typedef value_t*                                  value_ptr;

    /// type for head and tail position within block
    typedef ucl::aligned_atomic_type<value_ptr>       data_ptr;

    /// pointer to next data block
    typedef ucl::aligned_atomic_type<dataQ<value_t>*> next_ptr;

    /// number of elements in a block
    static constexpr nat_t BLK = DEFAULT_BLOCK_SIZE/2;

    data_ptr   hd;         ///< block internal head
    data_ptr   tl;         ///< block internal tail
    next_ptr   next;       ///< next block

    /// elements in the block
    /// \note in order to avoid interference from ctor and dtor calls
    ///       by classes with non-trivial constructor and deconstructor, the
    ///       block manages construction and deconstruction of its task objects.
    ///       Thus, it just reserves storage for BLK elements, and calls ctor
    ///       and dtor upon enqueue and dequeue (move based when possible).
    char       data[BLK*sizeof(value_t)];

    /// constructs data block with 0 elements
    dataQ()
    : hd(reinterpret_cast<value_ptr>(&data)),
      tl(reinterpret_cast<value_ptr>(&data)),
      next(nullptr)
    {}

    /// constructs a data block with 1 element
    /// \param el the first task
    explicit
    dataQ(value_t&& el)
    : hd(space_start()), tl(space_start() + 1), next(nullptr)
    {
      new (space_start()) value_t (std::move(el));
    }

    /// constructs a data block with 1 element
    /// \param el the first task
    explicit
    dataQ(const value_t& el)
    : hd(space_start()), tl(space_start() + 1), next(nullptr)
    {
      new (space_start()) value_t (el);
    }

    /// \private
    /// manage data space
    value_t* space_start()
    {
      return reinterpret_cast<value_t*>(&data);
    }

    /// \private
    /// manage data space
    const value_t* space_limit() const
    {
      return reinterpret_cast<const value_t*>(&data) + BLK;
    }

    /// tests if the data block can store more tasks
    bool has_space() const
    {
      // \mo relaxed, b/c tl is only modified by this thread
      return tl.val.load(std::memory_order_relaxed) < space_limit();
    }

    /// copy-enqueues a new element in this block
    /// \pre hasSpace() == true
    /// \pre only executed by the queue owner
    void enq(const value_t& el)
    {
      // \mo relaxed, b/c tl is only modified by this thread
      const value_ptr tail = tl.val.load(std::memory_order_relaxed);
      assert(tail < space_limit());

      new (tail) value_t (el); // in-place constructions

      // \mo release b/c data[tail] is published
      tl.val.store(tail+1, std::memory_order_release);
    }

    /// move-enqueues a new element in this block
    /// \pre hasSpace() == true
    /// \pre only executed by the queue owner
    void enq(value_t&& el)
    {
      // \mo relaxed, b/c tl is only modified by this thread
      const value_ptr tail = tl.val.load(std::memory_order_relaxed);
      assert(tail < space_limit());

      new (tail) value_t (std::move(el)); // move construction in-place

      // \mo release b/c data[tail] is published
      tl.val.store(tail+1, std::memory_order_release);
    }

    /// dequeues a task and moves it into res
    /// \result the block from where the element was dequeued;
    ///         nullptr if unsuccessful.
    dataQ<T>* deq(void* rawmem)
    {
      // \mo relaxed, b/c hd does not publish anything
      value_ptr head = hd.val.load(std::memory_order_relaxed);
      // \mo relaxed, b/c we are on the same thread
      value_ptr tail = tl.val.load(std::memory_order_relaxed);

      for (;;)
      {
        //~ if (head >= space_limit())
        if (head >= tail)
        {
          // \mo relaxed, intra thread dependency
          dataQ<T>* nxt = next.val.load(std::memory_order_relaxed);

          if (nxt == nullptr) return nullptr;

          return nxt->deq(rawmem);
        }

        //~ assert(head < space_limit());
        //~ if (head >= tail) return nullptr;
/*        
        bool succ = false;
        
        if (tail == space_limit())
        {
          head = hd.val.fetch_add(1, std::memory_order_relaxed);
          succ = head < space_limit();
        }
        else
        {
          succ = hd.val.compare_exchange_weak(head, head+1, std::memory_order_relaxed, std::memory_order_relaxed);
        }
        
        if (succ)
        {
          new (rawmem) value_t (std::move(*head));
          head->~value_t();
          return this;
        }
*/
        // \mo relaxed (succ) since the data is never rewritten
        // \mo relaxed (fail) b/c nothing happened
        if (hd.val.compare_exchange_weak(head, head+1, std::memory_order_relaxed, std::memory_order_relaxed))
        {
          assert(head < tail);

          // mv and deconstruct obj
          new (rawmem) value_t (std::move(*head));
          head->~value_t();
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
    /// \result success
    /// \pre    this is pinned, if a finegrain (i.e., publishAndScan) is used
    /// \post   this is unpinned
    /// \note   publishAndScan is supported, though finegrain memory management
    ///         does not bring anything to the table for this task implementation.
    template <class Alloc>
    bool deq_steal(void* rawmem, uint_fast8_t tries, Alloc alloc)
    {
      uint_fast8_t attempts = tries;

      // \mo relaxed, b/c hd does not publish anything
      value_ptr    head = hd.val.load(std::memory_order_relaxed);
      // \mo acquire, b/c we read from data[i], where i < tl
      value_ptr    tail = tl.val.load(std::memory_order_acquire);

      // give up under high contention
      while (attempts > 0)
      {
        //~ if (head == space_limit())
        if (head >= tail)
        {
          // \mo consume, b/c we want to see initialized memory
          dataQ<T>* nxt = alloc.template pin<std::memory_order_consume>(next.val);

          alloc.unpin(this, -1);
          if (nxt == nullptr) return false;

#if UCL_RUNTIME_DATA
          ++ucl::task_steal_skips;
#endif /* UCL_RUNTIME_DATA */
          return nxt->deq_steal(rawmem, tries, alloc);
          //~ return (nxt != nullptr && nxt->deq_steal(rawmem, tries, alloc));
        }

        //~ assert(head < space_limit());
        //~ if (head >= tail) return false;

        // \mo relaxed (succ) since the data is never rewritten
        // \mo relaxed (fail) b/c nothing happened
        if (hd.val.compare_exchange_weak(head, head+1, std::memory_order_relaxed, std::memory_order_relaxed))
        {
          assert(head+1 <= tail);

          new (rawmem) value_t(std::move(*head));
          head->~value_t();
          return true;
        }

        --attempts;
      }

      return false;
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
      static const nat_t                               SZPINWALL = 2;

    // private:
      /// The allocator/memory manager
      allocator_type                                   nodealloc;

      /// tail (enqueue-side) of the queue; only accessed by owner
      dataQ<T>*                                        tail; // only accessed from queue owner

      /// head (dequeue-side) of the queue; concurrently accessed by work stealers.
      ucl::aligned_atomic_type<dataQ<T>*, CACHELINESZ> head;

      /// number of tasks stolen from this thread
      nat_t                                            stolen;

      /// number of tasks enqueued by this thread
      nat_t                                            tasks_made;

      /// number of completed tasks (incl. tasks solved by other threads, but
      ///   not tasks that this thread stole).
      nat_t                                            tasks_done;

      /// number of stolen and completed tasks; updated by the work stealers.
      ucl::aligned_atomic_type<nat_t, CACHELINESZ>     steals;

      /// estimate of available tasks
      ucl::aligned_atomic_type<nat_t, CACHELINESZ>     numtasks;

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
      nat_t estimate_size() const
      {
        return numtasks.val.load(std::memory_order_relaxed);
      }

      /// allocates a new dataQ block
      dataQ<T>* new_block(T&& el)
      {
        allocator_type  alloc = get_allocator();
        dataQ<T>*       tmp   = alloc.allocate(1);

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
      bool deq(void* taskmem)
      {
        // \mo relaxed, b/c head is only updated by the owner
        dataQ<T>* curr = head.val.load(std::memory_order_relaxed);
        assert(curr);

        // attempt to dequeue from own queue
        dataQ<T>* actl = curr->deq(taskmem);

        // all blocks that have become empty can be freed
        if (curr != actl)
        {
#if UCL_RUNTIME_DATA     
          ucl_runtime_data before;
      
          ucl_rt_begin(before); 
#endif /* UCL_RUNTIME_DATA */     

          // \mo release to make content of new head available (needed?)
          allocator_type alloc = get_allocator();
          dataQ<T>*      limit = actl ? actl : tail;
          
          head.val.store(limit, std::memory_order_release);

          while (curr != limit)
          {
            dataQ<T>* tmp = curr;

            curr = curr->next.val.load(std::memory_order_relaxed);
            alloc.deallocate(tmp, 1);
          }

          alloc.release_memory_if_needed();

#if INFORMED_STEALING
          // update completed stolen tasks
          has_unfinished_tasks();
#endif /* INFORMED_STEALING */

          // update available task counter
          numtasks.val.store(tasks_made - tasks_done, std::memory_order_relaxed);

#if UCL_RUNTIME_DATA     
          ucl_rt_end(tk_deq_collect_data, before); 
#endif /* UCL_RUNTIME_DATA */     
        }

        return actl != nullptr;
      }

      /// attempts to steal a task from some other thread's queue
      /// \param  res   where the result will be stored
      /// \param  tries number of attempts before yielding to contention
      /// \result true if successful, false otherwise
      bool deq_steal(void* rawmem, uint_fast8_t tries)
      {
        // \mo consume, b/c we need to see memory initialized
        allocator_type  alloc = get_allocator();
        dataQ<T>*       curr  = alloc.template pin<std::memory_order_consume>(head.val);
        assert(curr);

        return curr->deq_steal(rawmem, tries, alloc);
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
      /// \param isMainThread The main thread will be responsible to
      ///                     reset all global state in the allocator.
      /// \details
      ///    DO NOT CALL when concurrent operations are ongoing!
      void qshutdown(bool isMainThread)
      {
        assert (head.val.load(std::memory_order_relaxed) == tail);
        
        get_allocator().deallocate(tail, 1);
        get_allocator().qrelease_memory();
        get_allocator().qdestruct(isMainThread);
      }
  };

  template <class T, class UclAlloc>
  bool flatQ<T, UclAlloc>::has_unfinished_tasks()
  {
    nat_t prev_stolen = stolen;

    stolen = steals.val.load(std::memory_order_relaxed);

    assert(tasks_made - tasks_done >= stolen - prev_stolen);
    tasks_done += (stolen - prev_stolen);

    return tasks_made > tasks_done;
  }


  template <class Q>
  struct active_tracker : Q
  {
    /// keeps track of current
    bool  active;

#if UCL_RUNTIME_DATA
    nat_t task_steals      = 0;
    nat_t task_steal_tries = 0;
    nat_t task_steal_skips = 0;
    
    ucl_runtime_data task_enq_data;
    ucl_runtime_data task_deq_data;
    ucl_runtime_data task_deq_steal_data;
    ucl_runtime_data task_deq_collect_data;
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
    typedef int_fast8_t                                  steal_try_cnt_t;

    /// maximum number of threads
    const nat_t                                          MAXTQ;

    /// the allocator/memory manager
    allocator_type                                       nodealloc;

    /// counts active threads; when count reaches 0 all task have been handled
    ucl::aligned_atomic_type<nat_t, CACHELINESZ>         active;

    /// task specific queue handles
    ucl::aligned_atomic_type<tasq*, CACHELINESZ>         taskq[256]; // \todo remove magic constant

    /// thread id within pool
    static thread_local nat_t                            idx;

    /// victim (thread id) where the last task has been stolen
    static thread_local nat_t                            last_victim;

    /// this task queue
    static thread_local tasq*                            tq_loc;

    /// last victim's task queue
    static thread_local tasq*                            tq_rem;

    /// initializes empty pool
    /// \param numthreads number of worker threads (and task queues)
    /// \param alloc      memory manager
    explicit
    perthread_pool(nat_t numthreads, const allocator_type& alloc = allocator_type())
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

      // set tasqs to null
      for (nat_t i = 1; i < MAXTQ; ++i)
      {
        // \mo relaxed, b/c the pool is constructed before we fork
        taskq[i].val.store(nullptr, std::memory_order_relaxed);
      }
    }

    /// initializes pool and enqueues first task
    /// \param numthreads number of worker threads (and task queues)
    /// \param work       first task
    /// \param alloc      memory manager
    explicit
    perthread_pool(nat_t numthreads, value_type work, const allocator_type& alloc = allocator_type())
    : perthread_pool(numthreads, alloc)
    {
      enq(std::move(work));
    }

    ~perthread_pool()
    {
#if UCL_RUNTIME_DATA
      // otherwise task pools are deleted at the end of the taskloop
      for (nat_t i = 0; i < MAXTQ; ++i)
      {
        delete taskq[i].val.load(std::memory_order_relaxed);
      }
#endif /* UCL_RUNTIME_DATA */
    }


    /// retrieve allocator
    allocator_type get_allocator()
    {
      return nodealloc;
    }

    /// returns the maximum number of worker threads for this pool
    nat_t num_threads() const { return MAXTQ; }

    /// initializes thread local storage of task pool
    /// \param parent pool id of parent thread
    /// \param self   pool id of this thread
    /// \note pool ids are generated by the task system
    void init(nat_t parent, nat_t self)
    {
      if (self == 0) return;

      assert(tq_loc == nullptr);

      // \todo \new
      tq_loc      = new tasq(nodealloc);
      idx         = self;
      last_victim = parent;

      taskq[self].val.store(tq_loc, std::memory_order_release);
    }

    /// copy-enqueues a new task
    void enq(const value_type& el)
    {
#if UCL_RUNTIME_DATA     
      ucl_runtime_data before;
        
      ucl_rt_begin(before); 
#endif /* UCL_RUNTIME_DATA */     
      
      tq_loc->enq(el);

#if UCL_RUNTIME_DATA     
      ucl_rt_end(tk_enq_data, before); 
#endif /* UCL_RUNTIME_DATA */     
    }

    /// move-enqueues a new task
    void enq(value_type&& el)
    {
#if UCL_RUNTIME_DATA     
      ucl_runtime_data before;
        
      ucl_rt_begin(before); 
#endif /* UCL_RUNTIME_DATA */     

      tq_loc->enq(std::move(el));

#if UCL_RUNTIME_DATA     
      ucl_rt_end(tk_enq_data, before); 
#endif /* UCL_RUNTIME_DATA */
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
      nat_t cntactv = became_inactive(tq_loc)
                          ? active.val.fetch_sub(1, std::memory_order_relaxed) - 1
                          : active.val.load(std::memory_order_relaxed)
                          ;

      return cntactv > 0;
    }

    /// dequeues a task and moves it into res
    /// \result  true if a task was dequeued, false otherwise
    bool deq(void* rawmem)
    {
#if UCL_RUNTIME_DATA     
      ucl_runtime_data before;
      
      ucl_rt_begin(before); 
#endif /* UCL_RUNTIME_DATA */     

      // try to dequeue from the local queue first
      if (tq_loc->deq(rawmem)) 
      { 
#if UCL_RUNTIME_DATA     
        ucl_rt_end(tk_deq_data, before); 
#endif /* UCL_RUNTIME_DATA */     
        
        return true;
      }

#if UCL_RUNTIME_DATA     
      ucl_rt_end(tk_deq_data, before); 
#endif /* UCL_RUNTIME_DATA */     

      // no more tasks in the local queue
      //   -> steal a task from somewhere else
#if UCL_RUNTIME_DATA     
      ucl_runtime_data before_steal;
      
      ucl_rt_begin(before_steal); 
#endif /* UCL_RUNTIME_DATA */     
      

#if INFORMED_STEALING
      // estimate average number of tasks from a sample
      //   and choose first victim from the sample with the most tasks.
      static const nat_t SAMPLESIZE = 4;

      nat_t            tasks_sum = 0;
      nat_t            tasks_avail[SAMPLESIZE];

      {
        nat_t          tasks_i = SAMPLESIZE;
        // nat_t          probe[SAMPLESIZE] = { last_victim, last_victim+1, idx+(MAXTQ-1), idx+1 };
        nat_t          probe[SAMPLESIZE] = { last_victim, last_victim, idx^2, idx^2 };
        nat_t          max     = 0;

        while (tasks_i)
        {
          --tasks_i;
          nat_t tid = probe[tasks_i];
          while (tid >= MAXTQ) { tid -= MAXTQ; }

          const tasq* victim = taskq[tid].val.load(std::memory_order_consume);

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

        allocator_type alloc = get_allocator();
        PinGuard       pinguard(alloc, SPMCBag::SZPINWALL);

        // work stealing
        nat_t          thrid   = last_victim;

#if INFORMED_STEALING
        nat_t          tasks_i = 0;
#endif /* INFORMED_STEALING */

        do
        {
          // \mo consume, b/c we follow the pointer
          tasq*         victim = taskq[thrid].val.load(std::memory_order_consume);

          if (victim)
          {
            int_fast8_t  tries = arch_model::num_tries(idx, thrid);

#if INFORMED_STEALING
            nat_t       avail = victim->estimate_size();

            if (avail * (SAMPLESIZE/2) > tasks_sum)
            {
              tries = std::min<int_fast8_t>(tries+4, 8);
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

            if (tries)
            {
              if (victim->deq_steal(rawmem, tries))
              {
                assert(idx != thrid);
#if UCL_RUNTIME_DATA
                ++task_steals;
                ucl_rt_end(tk_deq_steal_data, before_steal); 
#endif /* UCL_RUNTIME_DATA */

                tq_rem      = victim;
                last_victim = thrid;
                
                return true;
              }
            }

#if UCL_RUNTIME_DATA
            task_steal_tries += tries;
#endif /* UCL_RUNTIME_DATA */
          }

          ++thrid;
          if (thrid == MAXTQ) thrid = 0;
        } while (thrid != last_victim);
      }

      last_victim = 0;
#if UCL_RUNTIME_DATA
      ++task_steals;
      ucl_rt_end(tk_deq_steal_data, before_steal); 
#endif /* UCL_RUNTIME_DATA */
      return false;
    }

    /// releases all thread local memory
    void qshutdown(bool isMainThread)
    {
      tq_loc->qshutdown(isMainThread);
#if UCL_RUNTIME_DATA
      tq_loc->task_steals           = ucl::task_steals;
      tq_loc->task_steal_tries      = ucl::task_steal_tries;
      tq_loc->task_steal_skips      = ucl::task_steal_skips;
      tq_loc->task_enq_data         = ucl::tk_enq_data;
      tq_loc->task_deq_data         = ucl::tk_deq_data; 
      tq_loc->task_deq_steal_data   = ucl::tk_deq_steal_data;
      tq_loc->task_deq_collect_data = ucl::tk_deq_collect_data; 
#else
      delete tq_loc;
#endif /* UCL_RUNTIME_DATA */
    }
  };

  template <class Q>
  thread_local
  typename ucl::nat_t perthread_pool<Q>::idx = 0;

  template <class Q>
  thread_local
  typename ucl::nat_t perthread_pool<Q>::last_victim = 0;

  template <class Q>
  thread_local
  typename perthread_pool<Q>::tasq* perthread_pool<Q>::tq_loc = nullptr;

  template <class Q>
  thread_local
  typename perthread_pool<Q>::tasq* perthread_pool<Q>::tq_rem = nullptr;

  /// \private
  /// \brief cleanup after all threads finish
  template <class P>
  void threadexit(std::atomic<nat_t>& signal, std::atomic<nat_t>& coordinator, P& taskpool)
  {
    // read current barrier value
    const nat_t eofwork_barrier_value = coordinator.load(std::memory_order_relaxed);
    
    // notify parent about finished computation
    signal.store(1, std::memory_order_release);

    // wait until all threads are confirmed to be in quiesence state
    while (eofwork_barrier_value == coordinator.load(std::memory_order_acquire));

    // release resources
    taskpool.qshutdown(false);
  }

  /// \private
  /// \brief  main task loop
  /// \tparam R result type
  /// \tparam P pool type
  /// \tparam F function type
  template <class R, class P, class F>
  static
  void taskloop(nat_t parent, nat_t thrnum, R* res, P* tasks, F fun)
  {
    typedef typename P::value_type                                                     task_type;
    typedef typename std::aligned_storage<sizeof(task_type), alignof(task_type)>::type aligned_mem_type;

#if !defined __CYGWIN__ && !defined OS_WIN32 && !defined __OpenBSD__ && !defined __sun
    arch_model::bind_to_core(pthread_self(), thrnum);
#endif

    assert(res && tasks);
    assert(thrnum < tasks->num_threads());

    // initialize task-queues for each thread
    tasks->init(parent, thrnum);

    // create task storage
    aligned_mem_type work_space;
    task_type* const work = reinterpret_cast<task_type*>(&work_space);

    // initialize local reduction variable
    R                val = R{};

    // run until no more tasks can be found
    do
    {
      while (tasks->deq(work))
      {
        tasks->work_started();
        val += fun(*tasks, std::move(*work));
        work->~task_type();
        tasks->work_completed();
      }
    } while (tasks->work_available());

    *res = val;
  }

  /// \private
  /// \brief calls main task loop and manages resource upon termination
  template <class R, class P, class F>
  static
  void taskloop_term( std::atomic<nat_t>* signal, std::atomic<nat_t>* coordinator,
                      nat_t parent, nat_t thrnum,
                      R* res, P* tasks, F fun
                    )
  {
    assert(signal && coordinator);

    taskloop<R, P, F>(parent, thrnum, res, tasks, fun);

    threadexit(*signal, *coordinator, *tasks);
  }

  /// \private
  /// \brief  spawns n threads sequentially
  /// \tparam R result type
  /// \tparam P pool type
  /// \tparam F function type
  template <class R, class P, class F>
  static
  void spawn_threads_single( std::atomic<nat_t>* coordinator,
                             nat_t parent, nat_t thrnumlo, nat_t thrnumhi,
                             R* res, P* taskpool, F worker
                           )
  {
    assert(thrnumhi > 0);
    const nat_t last = thrnumhi-1;

    if (thrnumlo == last)
    {
      taskloop(parent, thrnumlo, res, taskpool, worker);
      return;
    }

    R                  sub;
    std::atomic<nat_t> signal(0);
    ucl::thread        t(taskloop_term<R,P,F>, &signal, coordinator, thrnumlo, last, &sub, taskpool, worker);

    t.detach();
    spawn_threads_single<R>(coordinator, thrnumlo, thrnumlo, last, res, taskpool, worker);
    while (!signal.load(std::memory_order_acquire));

    *res += sub;
  }

  template <class R, class P, class F>
  static
  void spawn_threads_term( std::atomic<nat_t>* parent_signal, std::atomic<nat_t>* coordinator,
                           nat_t parent, nat_t thrnumlo, nat_t thrnumhi,
                           R* res, P* taskpool, F worker
                         );


  /// \private
  /// \brief  spawns threads; if the number exceeds a preset threshold
  ///         the thread spawns another spawning thread before spawning
  ///         half the threads.
  /// \tparam R result type
  /// \tparam P pool type
  /// \tparam F function type
  template <class R, class P, class F>
  static
  void spawn_threads( std::atomic<nat_t>* coordinator,
                      nat_t parent, nat_t thrnumlo, nat_t thrnumhi,
                      R* res, P* taskpool, F worker
                    )
  {
    static const nat_t SPAWN_THREASHOLD = 8;

    const nat_t numthreads = thrnumhi - thrnumlo;

    if (numthreads <= SPAWN_THREASHOLD)
    {
      spawn_threads_single(coordinator, parent, thrnumlo, thrnumhi, res, taskpool, worker);
      return;
    }

    const nat_t        thrnummid = thrnumlo+numthreads/2;
    R                  sub;
    std::atomic<nat_t> signal(0);
    ucl::thread        spawner(spawn_threads_term<R,P,F>, &signal, coordinator, thrnumlo, thrnummid, thrnumhi, &sub, taskpool, worker);

    spawner.detach();
    spawn_threads<R>(coordinator, thrnumlo, thrnumlo, thrnummid, res, taskpool, worker);
    while (!signal.load(std::memory_order_acquire));

    *res += sub;
  }

  template <class R, class P, class F>
  static
  void spawn_threads_term( std::atomic<nat_t>* parent_signal, std::atomic<nat_t>* coordinator,
                           nat_t parent, nat_t thrnumlo, nat_t thrnumhi,
                           R* res, P* taskpool, F worker
                         )
  {
    spawn_threads<R>(coordinator, parent, thrnumlo, thrnumhi, res, taskpool, worker);

    threadexit(*parent_signal, *coordinator, *taskpool);
  }


/*
  /// \brief  spawns a spawning thread and continues to do work
  /// \tparam R result type
  /// \tparam P pool type
  /// \tparam F function type
  template <class R, class P, class F>
  static
  void spawn_spawner(nat_t thrnumlo, nat_t thrnumhi, R* res, P* taskpool, F worker)
  {
    static const nat_t MIN_THREADS = 2;

    const nat_t numthreads = thrnumhi - thrnumlo;

    if (numthreads < MIN_THREADS)
    {
      spawn_threads_single(thrnumlo, thrnumhi, res, taskpool, worker);
      return;
    }

    const nat_t thrnummid = thrnumlo+1;
    R            sub;
    ucl::thread  spawner(spawn_threads<R,P,F>, thrnummid, thrnumhi, &sub, taskpool, worker);

    taskloop(thrnumlo, res, taskpool, worker);
    spawner.join();

    *res += sub;
  }
*/

#if UCL_RUNTIME_DATA
  template <class T>
  static inline
  void log_runtime_data(std::ostream& os, ucl::perthread_pool<T>& pool)
  {
    typedef typename ucl::perthread_pool<T>::tasq tasq;

    nat_t total_stolen           = 0;
    nat_t total_tasks            = 0;
    nat_t total_task_steals      = 0;
    nat_t total_task_steal_tries = 0;
    nat_t total_task_steal_skips = 0;

    for (nat_t i = 0; i < pool.MAXTQ; ++i)
    {
      tasq& tsk = *pool.taskq[i].val;

      total_tasks            += tsk.tasks_done;
      total_stolen           += tsk.stolen;
      total_task_steals      += tsk.task_steals;
      total_task_steal_tries += tsk.task_steal_tries;
      total_task_steal_skips += tsk.task_steal_skips;

      os << "t " << i
         << "  " << tsk.tasks_done
         << " (" << tsk.stolen << ")"
         << " [" << tsk.task_steals << "/" << tsk.task_steal_tries << "/" << tsk.task_steal_skips << "]"
         << " <" << tsk.task_enq_data << "/" << tsk.task_deq_data 
         << "/"  << tsk.task_deq_collect_data << "/" << tsk.task_deq_steal_data
         << ">"  
         << std::endl;
    }

    os << "**** totals: " << total_tasks
       << " (" << total_stolen << " " << ((total_stolen * 100.0) / (total_tasks))
       << "%/" << total_task_steals
       << "/" << total_task_steal_tries << "/" << total_task_steal_skips
       << ")"
       << std::endl;

    log_epoch_telemetry(std::cerr, pool.get_allocator());
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
       -> decltype( fun(pool, std::declval<typename TaskPool::value_type>()) )
  {
    using R = decltype( fun(pool, std::declval<typename TaskPool::value_type>()) );

    static std::atomic<nat_t> coordinator(0);
    
    const nat_t numthreads = pool.num_threads();
    R           res;

    arch_model::set_threadinfo_info(numthreads);

    spawn_threads(&coordinator, 0, 0, numthreads, &res, &pool, fun);

    coordinator.fetch_add(1, std::memory_order_release);
    pool.qshutdown(true);

#if UCL_RUNTIME_DATA
    log_runtime_data(std::cerr, pool);
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
  auto execute_tasks(nat_t numthreads, F fun, T&& task)
       -> decltype(execute_pool(fun, *new fifo_queue_pool<T>(numthreads, std::forward<T>(task))))
  {
    fifo_queue_pool<T> pool(numthreads, std::forward<T>(task));

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
  auto execute_tasks(nat_t numthreads, F fun, T&& task)
       -> decltype(execute_pool(fun, *new perthread_pool<T,A>(numthreads, std::move(task))))
  {
    perthread_pool<T,A> taskpool(numthreads, std::move(task));

    return execute_pool(fun, taskpool);
  }
#endif

} // end namespace ucl

#endif /* _UCL_TASK_HPP */
