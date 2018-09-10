#ifndef _TASKS_HPP
#define _TASKS_HPP

#include <algorithm>
#include <thread>

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
  static const size_t BLKSZ                = 1024;

  template <class T>
  struct dataQ
  {
    typedef uab::aligned_atomic_type<size_t, CACHELINESZ>    counter_t;
    typedef uab::aligned_atomic_type<dataQ<T>*, CACHELINESZ> next_t;
    typedef T                                                value_t;

    static const size_t BLK = BLKSZ;

    counter_t  hd;
    counter_t  tl;
    next_t     next;
    size_t     dummy;
    value_t    data[BLK];

    dataQ()
    : hd(0), tl(0), next(nullptr)
    {}

    explicit
    dataQ(T el)
    : hd(0), tl(1), next(nullptr)
    {
      data[0] = el;
    }

    bool has_space()
    {
      // \mo relaxed, b/c tl is only modified by this thread
      return tl.val.load(std::memory_order_relaxed) < BLK;
    }

    void enq(const T& el)
    {
      // \mo relaxed, b/c tl is only modified by this thread
      const size_t tail = tl.val.load(std::memory_order_relaxed);

      assert(tail < BLK);
      data[tail] = el;

      // \mo release b/c data[tail] is published
      tl.val.store(tail+1, std::memory_order_release);
    }

    dataQ<T>* deq(T& res)
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

        res = data[head];
        // \mo release to prevent reordering with previous read of data[head]
        //   \todo since the data is not rewritten, relaxed might do
        // \mo relaxed b/c nothing happened
        if (hd.val.compare_exchange_strong(head, head+1, std::memory_order_relaxed, std::memory_order_relaxed))
        {
          assert(head+1 <= tail);
          return this;
        }
      }
    }

    template <class Alloc>
    dataQ<T>* deq_steal(T& res, size_t tries, Alloc alloc)
    {
      size_t attempts = tries;

      // \mo relaxed, b/c hd does not publish anything
      size_t head = hd.val.load(std::memory_order_relaxed);
      // \mo acquire, b/c we read from data[i], where i < tl
      size_t tail = tl.val.load(std::memory_order_acquire);

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

        res = data[head];
        // \mo release to prevent reordering with previous read of data[head]
        //   \todo since the data is not rewritten, relaxed might do
        // \mo relaxed b/c nothing happened
        if (hd.val.compare_exchange_strong(head, head+1, std::memory_order_relaxed, std::memory_order_relaxed))
        {
          assert(head+1 <= tail);
          return this;
        }

        --attempts;
      }

      return nullptr;
    }
  };


  template <class T, class Alloc>
  struct flatQ
  {
    typedef Alloc                                                        allocator_type;
    typedef std::allocator_traits<Alloc>                                 orig_alloc_traits;
    typedef typename orig_alloc_traits::template rebind_alloc<dataQ<T> > node_alloc_type;

    static const size_t                              SZPINWALL = 2;

    node_alloc_type                                  nodealloc;
    bool                                             active;
    size_t                                           stolen;
    size_t                                           tasks_made;
    size_t                                           tasks_done;
    dataQ<T>*                                        tail; // only accessed from queue owner

    uab::aligned_atomic_type<dataQ<T>*, CACHELINESZ> head;
    uab::aligned_atomic_type<size_t, CACHELINESZ>    numtasks;
    uab::aligned_atomic_type<size_t, CACHELINESZ>    steals;

    explicit
    flatQ(allocator_type alloc = allocator_type())
    : nodealloc(alloc), active(true), stolen(0), tasks_made(0), tasks_done(0),
      tail(new dataQ<T>()), head(tail), numtasks(0), steals(0)
    {}

    node_alloc_type get_allocator()
    {
      return nodealloc;
    }

    dataQ<T>* new_tile(T el)
    {
      node_alloc_type alloc = get_allocator();
      dataQ<T>*       tmp = alloc.allocate(1);

      alloc.construct(tmp, el);
      return tmp;
    }

    void enq(const T& el)
    {
      ++tasks_made;

      if (tail->has_space())
      {
        tail->enq(el);
        return;
      }

      dataQ<T>* tmp = new_tile(el);

      numtasks.val.store(tasks_made - tasks_done, std::memory_order_relaxed);

      // \mo release, b/c new_tile was initialized
      tail->next.val.store(tmp, std::memory_order_release);
      tail = tmp;
    }

    void task_completed()
    {
      ++tasks_done;
    }

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

    bool deq(T& res)
    {
      typedef typename node_alloc_type::pinguard PinGuard;

      // \mo relaxed, b/c head is only updated by the owner
      dataQ<T>*               curr = head.val.load(std::memory_order_relaxed);
      assert(curr);

      dataQ<T>* actl = curr->deq(res);

      if (actl == nullptr) return false;

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

        // \mo release to make content of new head available
        head.val.store(actl, std::memory_order_release);
        numtasks.val.store(tasks_made - tasks_done, std::memory_order_relaxed);
      }

      return true;
    }

    bool deq_steal(T& res, size_t tries)
    {
      typedef typename node_alloc_type::pinguard PinGuard;

      node_alloc_type alloc = get_allocator();
      PinGuard        pinguard(alloc, SZPINWALL);

      // \mo consume, b/c we need to see memory initialized
      dataQ<T>*       curr = alloc.template pin<std::memory_order_consume>(head.val);
      assert(curr);

      return curr->deq_steal(res, tries, alloc) != nullptr;
    }

    void qrelease_memory()
    {
      get_allocator().qrelease_memory();
    }
  };


  template <class T, class _Alloc = lockfree::epoch_manager<T, std::allocator>>
  struct pool
  {
    typedef T                                                       task_type;
    typedef _Alloc                                                  allocator_type;
    typedef flatQ<T, _Alloc>                                        tasq;

    typedef std::allocator_traits<_Alloc>                           orig_alloc_traits;
    typedef typename orig_alloc_traits::template rebind_alloc<tasq> node_alloc_type;

    const size_t                                  MAXTQ;

    node_alloc_type                               nodealloc;

    uab::aligned_atomic_type<size_t, CACHELINESZ> active;
    uab::aligned_atomic_type<tasq*, CACHELINESZ>  taskq[256]; // \todo remove magic constant

    static thread_local size_t                    idx;
    static thread_local size_t                    last_victim;
    static thread_local tasq*                     tq_loc;
    static thread_local tasq*                     tq_rem;

    explicit
    pool(size_t numthreads, T work, const allocator_type& alloc = allocator_type())
    : MAXTQ(numthreads), nodealloc(alloc), active(numthreads)
    {
      assert(tq_loc == nullptr);
      assert(MAXTQ <= 256);

      // initialize main thread
      tq_loc      = new tasq(nodealloc);
      idx         = 0;
      last_victim = 0;
      taskq[0].val.store(tq_loc, std::memory_order_relaxed);

      enq(work);

      // set tasqs to null
      for (size_t i = 1; i < MAXTQ; ++i)
      {
        // \mo relaxed, b/c the pool is constructed before we fork
        taskq[i].val.store(nullptr, std::memory_order_relaxed);
      }
    }

    void init(size_t parent, size_t self)
    {
      if (self == 0) return;

      tq_loc      = new tasq(nodealloc);
      idx         = self;
      last_victim = parent;

      taskq[self].val.store(tq_loc, std::memory_order_release);
    }

    void enq(T el)
    {
      tq_loc->enq(el);
    }

    void work_started()
    {
      if (!tq_loc->active)
      {
        active.val.fetch_add(1, std::memory_order_relaxed);
        tq_loc->active = true;
      }
    }

    void work_completed()
    {
      if (tq_rem == nullptr)
      {
        tq_loc->task_completed();
        return;
      }

      tq_rem->steals.val.fetch_add(1, std::memory_order_relaxed);
      tq_rem = nullptr;
    }

    void check_active()
    {
      if (tq_loc->became_inactive())
      {
        active.val.fetch_sub(1, std::memory_order_relaxed);
      }
    }

    bool has_work()
    {
      return active.val.load(std::memory_order_relaxed) > 0;
    }

    bool deq(T& res)
    {
      static const size_t SAMPLESIZE = 4;

      bool succ = tq_loc->deq(res);

      if (succ) return succ;

      // estimate average number of tasks by using rolling average
      size_t             tasks_avail[SAMPLESIZE];
      size_t             tasks_sum = 0;
      size_t             thrid     = last_victim;

      {
        size_t           tasks_i = 0;
        size_t           tid     = thrid;
        size_t           max     = 0;

        while (tasks_i < SAMPLESIZE)
        {
          tasq*          victim = taskq[tid].val.load(std::memory_order_consume);

          if (victim)
          {
            tasks_avail[tasks_i] = victim->numtasks.val.load(std::memory_order_relaxed);
            tasks_sum += tasks_avail[tasks_i];

            if (tasks_avail[tasks_i] > max)
            {
              max = tasks_avail[tasks_i];
              last_victim = thrid = tid;
            }
          }

          ++tasks_i;
          tid += 3;
          while (tid >= MAXTQ) { tid -= MAXTQ; }
        }
      }


      {
        // work stealing
        size_t             tasks_i = 0;

        do
        {
          // \mo consume, b/c we follow the pointer
          tasq* victim = taskq[thrid].val.load(std::memory_order_consume);

          if (victim)
          {
            size_t tries = arch_model::num_tries(idx, thrid);
            size_t avail = victim->numtasks.val.load(std::memory_order_relaxed);

            if (avail * (SAMPLESIZE/2) > tasks_sum)
            {
              tries += 4;

              tries = std::min(tries, size_t(8));
            }
            else if (avail * (SAMPLESIZE*4) < tasks_sum)
            {
              // tries -=2;
              tries = 0;
            }

            tasks_sum -= tasks_avail[tasks_i];
            tasks_sum += (tasks_avail[tasks_i] = avail);
            ++tasks_i;
            if (tasks_i == SAMPLESIZE) tasks_i = 0;

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

    void qrelease_memory()
    {
      tq_loc->qrelease_memory();
    }
  };

  template <class W, class A>
  thread_local
  size_t pool<W,A>::idx = 0;

  template <class W, class A>
  thread_local
  size_t pool<W,A>::last_victim = 0;

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
    size_t core_init; // core where thread began
    size_t core_last;  // core where thread ended
  };

  static threadstat work[NUMTHREADS];
#endif /* PRINT_STATS */

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

  /// auxiliary class for non-blaze task systems
  template <class T>
  struct reducer
  {
    uab::aligned_type<T, CACHELINESZ>   myres[256];

    static std::atomic<int>             ctr;
    static thread_local size_t          val;
  };

  template <class T>
  std::atomic<int> reducer<T>::ctr(0);

  template <class T>
  thread_local
  size_t reducer<T>::val(ctr.fetch_add(1, std::memory_order_relaxed));

} // end namespace uab

#endif /* _TASKS_HPP */
