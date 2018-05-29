#ifndef _TASKS_HPP
#define _TASKS_HPP

#if PRINT_STATS
#include <sched.h>
#include <pthread.h>
#endif /* PRINT_STATS */

#include "atomicutil.hpp"
#include "pmemory.hpp"

// INTEL_DUAL_SOCKET POWER8_DUAL_SOCKET INTEL_PHI_1

#define INTEL_DUAL_SOCKET 0

namespace uab
{
  static const size_t BLKSZ           = 1024;

#if defined(INTEL_DUAL_SOCKET)
  typedef uab::aligned_type<size_t, CACHELINESZ> core_t;

  static const size_t NUMCORES = 40;

  static inline
  void bind_to_core(pthread_t thr, size_t num, size_t numthreads)
  {
    assert((num < NUMCORES) && (numthreads <= NUMCORES));

    // core 0,2,4,..,18,20,22,..,38: 0            <= num < numthreads/2 -> (2*num % numthreads)
    // core 1,3,5,..,19,21,23,..,39: numthreads/2 <= num < numthreads   -> (2*(num - (numthreads/2) ~1) % numthreads) | 1
    if (num < numthreads/2)
    {
      num = 2 * num;
    }
    else
    {
      num = (2 * (num - numthreads/2)) | 1;
    }

    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(num, &cpuset);

    /* int rc = */ pthread_setaffinity_np(thr, sizeof(cpu_set_t), &cpuset);
  }


  static inline
  size_t num_tries(size_t one, size_t two)
  {
    // steal from primary threads on each core
    if (two == 0) return 6;

    // same socket
    if (((one ^ two) & 1) == 0) return 4;

    // NUMA
    return 1;
  }


#elif defined(INTEL_PHI_1)

  static const size_t CORES   = 4;
  static const size_t CPUS    = 61;
  static const size_t THREADS = CORES * CPUS;

  static inline
  void bind_to_core(pthread_t thr, size_t num, size_t /*numthreads*/)
  {
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(num, &cpuset);

    /* int rc = */ pthread_setaffinity_np(thr, sizeof(cpu_set_t), &cpuset);
  }

  static inline
  size_t num_tries(size_t lhs, size_t rhs)
  {
    // same CPU
    if ((lhs & ~3) == (rhs & ~3) && ((lhs & 3) != 0))
    {
      return 4;
    }

    size_t dist = std::abs(int(lhs) - int(rhs));

    dist = std::min(dist, THREADS-dist);

    size_t res;

    if (dist <= 4 * CORES)       res = 2;
    else if (rhs == 0)           res = 1;
    else                         res = 0;

    if ((lhs & 3) == 0) res += 1;

    return res;
  }

#elif defined(POWER8_DUAL_SOCKET)

  static const size_t CORE_MULTIPLIER = 8; // no hyperthreading

  int bind_to_core(pthread_t thr, size_t num, size_t /*numthreads*/)
  {
    static_assert((CORE_MULTIPLIER & (CORE_MULTIPLIER-1)) == 0, "assumed power of 2");

    //  0 -> 0,  1 ->  8, 2 -> 16, .., 19 -> 152
    // 20 -> 1, 21 ->  9,
    // 40 -> 2, 41 -> 10,

    const size_t thread_on_core = num / 20;
    const size_t core_num       = num % 20;

    num = core_num * CORE_MULTIPLIER + thread_on_core;

    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(num, &cpuset);

    int rc = pthread_setaffinity_np(thr, sizeof(cpu_set_t), &cpuset);

    assert(rc == 0);
    return rc;
  }

  size_t num_tries(size_t thisthr, size_t thatthr)
  {
    if (thatthr == 0) return 8;

    size_t thiscore = thisthr / CORE_MULTIPLIER;
    size_t thatcore = thatthr / CORE_MULTIPLIER;

    if (thiscore == thatcore) return 8;

    if ((thiscore < 10) == (thatcore < 10)) return 4;

    return 1;
  }

#else /* DEFAULT */

  void bind_to_core(pthread_t, size_t, size_t) {}

  size_t num_tries(size_t, size_t two)
  {
    if (two == 0) return 6;

    return 2;
  }

#endif

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

    void enq(T el)
    {
      // \mo relaxed, b/c tl is only modified by this thread
      const size_t tail = tl.val.load(std::memory_order_relaxed);

      assert(tail < BLK);
      data[tail] = el;

      // \mo release b/c data[tail] is published
      tl.val.store(tail+1, std::memory_order_release);
    }

    std::pair<T, dataQ<T>*> deq()
    {
      T      res;
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

          if (nxt == nullptr)
          {
            return std::make_pair(T(), nullptr);
          }

          return nxt->deq();
        }

        assert(head < BLK);
        if (head >= tail)
        {
          return std::make_pair(T(), nullptr);
        }

        res = data[head];
        // \mo release to prevent reordering with previous read of data[head]
        //   \todo since the data is not rewritten, relaxed might do
        // \mo relaxed b/c nothing happened
        if (hd.val.compare_exchange_strong(head, head+1, std::memory_order_relaxed, std::memory_order_relaxed))
        {
          assert(head+1 <= tail);
          return std::make_pair(res, this);
        }
      }
    }

    template <class Alloc>
    std::pair<T, dataQ<T>*> deq_steal(size_t tries, Alloc alloc)
    {
      T      res;

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

          if (nxt == nullptr)
          {
            return std::make_pair(T(), nullptr);
          }

          return nxt->deq_steal(tries, alloc);
        }

        assert(head < BLK);
        if (head >= tail)
        {
          return std::make_pair(T(), nullptr);
        }

        res = data[head];
        // \mo release to prevent reordering with previous read of data[head]
        //   \todo since the data is not rewritten, relaxed might do
        // \mo relaxed b/c nothing happened
        if (hd.val.compare_exchange_strong(head, head+1, std::memory_order_relaxed, std::memory_order_relaxed))
        {
          assert(head+1 <= tail);
          return std::make_pair(res, this);
        }

        --attempts;
      }

      return std::make_pair(T(), nullptr);
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
    size_t                                           created;
    dataQ<T>*                                        tail; // only accessed from queue owner
    //~ dataQ<T>*                                        tile; // only accessed from queue owner

    uab::aligned_atomic_type<dataQ<T>*, CACHELINESZ> head;
    uab::aligned_atomic_type<size_t, CACHELINESZ>    steals;

    explicit
    flatQ(allocator_type alloc = allocator_type())
    : nodealloc(alloc), active(true), stolen(0), tasks_made(0), tasks_done(0),
      tail(new dataQ<T>()), head(tail), steals(0)
    {}

    node_alloc_type get_allocator()
    {
      return nodealloc;
    }

    dataQ<T>* new_tile(T el)
    {
      dataQ<T>* tmp = nullptr; // tile;

      if (tmp == nullptr)
      {
        node_alloc_type alloc = get_allocator();

        tmp = alloc.allocate(1);
        alloc.construct(tmp, el);
      }
      else
      {
        assert(false);
        tmp->~dataQ<T>();
        new (tmp) dataQ<T>(el);
        //~ tile = nullptr;
      }

      return tmp;
    }

    void enq(T el)
    {
      ++tasks_made;

      if (tail->has_space())
      {
        tail->enq(el);
        return;
      }

      dataQ<T>* tmp = new_tile(el);

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

    std::pair<T, bool> deq()
    {
      typedef typename node_alloc_type::pinguard PinGuard;

      // \mo relaxed, b/c head is only updated by the owner
      dataQ<T>*               curr = head.val.load(std::memory_order_relaxed);
      assert(curr);

      std::pair<T, dataQ<T>*> actl = curr->deq();

      if (actl.second == nullptr)
      {
        return std::make_pair(actl.first, false);
      }

      if (curr != actl.second)
      {
        node_alloc_type alloc = get_allocator();
        PinGuard        pinguard(alloc);

        do
        {
          dataQ<T>* tmp = curr;

          curr = curr->next.val.load(std::memory_order_relaxed);
          alloc.deallocate(tmp, 1);
        } while (curr != actl.second);

        // \mo release to make content of new head available
        head.val.store(actl.second, std::memory_order_release);
      }

      return std::make_pair(actl.first, true);
    }

    std::pair<T, bool> deq_steal(size_t tries)
    {
      typedef typename node_alloc_type::pinguard PinGuard;

      node_alloc_type alloc = get_allocator();
      PinGuard        pinguard(alloc, SZPINWALL);

      // \mo consume, b/c we need to see memory initialized
      dataQ<T>*               curr = alloc.template pin<std::memory_order_consume>(head.val);
      assert(curr);

      std::pair<T, dataQ<T>*> actl = curr->deq_steal(tries, alloc);

      return std::make_pair(actl.first, actl.second != nullptr);
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

      for (size_t i = 0; i < MAXTQ; ++i)
      {
        // \mo relaxed, b/c the pool is constructed before we fork
        taskq[i].val.store(nullptr, std::memory_order_relaxed);
      }

      assert(tq_loc == nullptr);
      tq_loc = new tasq();
      enq(work);
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

    std::pair<T, bool> deq()
    {
      std::pair<T, bool> res = tq_loc->deq();
      if (res.second) return res;

      // work stealing
      size_t             thrid  = last_victim;

      do
      {
        // \mo consume, b/c we follow the pointer
        tasq* victim = taskq[thrid].val.load(std::memory_order_consume);

        if (victim)
        {
          res = victim->deq_steal(num_tries(idx, thrid));

          if (res.second)
          {
            tq_rem      = victim;
            last_victim = thrid;
            return res;
          }
        }

        ++thrid;
        if (thrid == MAXTQ) thrid = 0;
      } while (thrid != last_victim);

      last_victim = 0;
      assert(!res.second);
      return res;
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

    bind_to_core(pthread_self(), thrnum, tasks->MAXTQ);

#if PRINT_STATS
    work[thrnum].core_init = sched_getcpu();
#endif /* PRINT_STATS */

    assert(thrnum < tasks->MAXTQ);

    // initialize task-queues for each thread
    if (P::tq_loc == nullptr)
      P::tq_loc = new typename P::tasq();

    P::idx = thrnum;
    P::last_victim = parent;

    tasks->taskq[thrnum].val.store(P::tq_loc, std::memory_order_release);

    // initialize local reduction variable
    R      val = R();

    // run until no more tasks can be found
    do
    {
      std::pair<task_type, bool> work = tasks->deq();

      while (work.second)
      {
        tasks->work_started();
        val  += fun(*tasks, work.first);
        tasks->work_completed();

        work  = tasks->deq();
      }

      tasks->check_active();
    } while (tasks->has_work());

    //~ std::pair<T, bool> test = pool<T>::tq_loc->deq();
    //~ assert(!test.second);

    *res = val;

#if PRINT_STATS
    work[thrnum].num       = P::tq_loc->tasks_done;
    work[thrnum].core_last = sched_getcpu();
#endif
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
