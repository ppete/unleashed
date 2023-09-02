
#ifndef _UCL_ARCHMODEL_H
#define _UCL_ARCHMODEL_H 1

#include <cassert>

#include "ucl/unused.hpp"
#include "ucl/thread.hpp"

//
// provides implementation of hardware models

namespace ucl
{
  size_t thread_mult; // == NUM_THREADS_PER_CORE if hyperthreading is not needed
  size_t thread_mask; // used to identify sibling hyperthreads

  struct generic_arch
  {
    static
    void set_threadinfo_info(size_t) {}

    static
    void bind_to_core(pthread_t, size_t) {}

    static
    size_t num_tries(size_t, size_t)
    {
      return 4;
    }

    size_t max_concurrency()
    {
      size_t hw = ucl::thread::hardware_concurrency();

      return hw ? hw:168;  /* artificial limit */
    }
  };

#if !defined __CYGWIN__ && !defined OS_WIN32 && !defined __OpenBSD__ && !defined __sun && !defined __APPLE__ && !defined _AIX
  template <size_t NSOCKETS, size_t NCORES, size_t NTHREADS>
  struct ascending
  {
    enum
    {
      NUM_SOCKETS            = NSOCKETS,
      NUM_CORES_PER_SOCKET   = NCORES,
      NUM_THREADS_PER_CORE   = NTHREADS,
      NUM_THREADS_PER_SOCKET = NUM_CORES_PER_SOCKET * NTHREADS,
      NUM_CORES              = NUM_SOCKETS * NUM_CORES_PER_SOCKET
    };

    static
    void set_threadinfo_info(size_t numthreads)
    {
      assert(numthreads <= NUM_CORES), ucl::unused(numthreads);
    }

    static
    void bind_to_core(pthread_t thr, size_t cpu)
    {
      cpu = cpu % NUM_CORES;

      cpu_set_t cpuset;
      CPU_ZERO(&cpuset);
      CPU_SET(cpu, &cpuset);
      CPU_SET(cpu+NUM_CORES, &cpuset);

      /* int rc = */ pthread_setaffinity_np(thr, sizeof(cpu_set_t), &cpuset);
    }

    static
    size_t num_tries(size_t thisthr, size_t thatthr)
    {
      const size_t thiscore = thisthr % NUM_CORES;
      const size_t thatcore = thatthr % NUM_CORES;

      if (thiscore == thatcore) return 6;

      const size_t thissock = thiscore / NUM_CORES_PER_SOCKET;
      const size_t thatsock = thatcore / NUM_CORES_PER_SOCKET;

      if (thissock == thatsock) return 4;

      return 3;
    }

    static constexpr
    size_t max_concurrency()
    {
      return NUM_SOCKETS * NUM_CORES_PER_SOCKET * NUM_THREADS_PER_CORE;
    }
  };

  template <size_t NSOCKETS, size_t NCORES, size_t NTHREADS>
  struct intel_arch
  {
    enum
    {
      NUM_SOCKETS            = NSOCKETS,
      NUM_CORES_PER_SOCKET   = NCORES,
      NUM_THREADS_PER_CORE   = NTHREADS,
      NUM_THREADS_PER_SOCKET = NUM_CORES_PER_SOCKET * NUM_THREADS_PER_CORE,
      NUM_CORES              = NUM_SOCKETS * NUM_CORES_PER_SOCKET
    };

    static
    void set_threadinfo_info(size_t numthreads)
    {
      assert(numthreads <= NUM_CORES);

      thread_mult = std::min<size_t>(NUM_CORES / numthreads, NUM_THREADS_PER_CORE);

      thread_mask = NUM_THREADS_PER_CORE / thread_mult;
      assert(thread_mask && (thread_mask & (thread_mask-1)) == 0);

      thread_mask = ~(thread_mask-1);
    }

    static
    void bind_to_core(pthread_t thr, size_t cpu)
    {
      assert(cpu < NUM_CORES);

      // core 0,2,4,..,18,20,22,..,38: 0            <= num < numthreads/2 -> (2*num % numthreads)
      // core 1,3,5,..,19,21,23,..,39: numthreads/2 <= num < numthreads   -> (2*(num - (numthreads/2) ~1) % numthreads) | 1

      const size_t socket_num = size_t(cpu >= (NUM_CORES_PER_SOCKET * NUM_THREADS_PER_CORE));

      cpu = ((cpu - (NUM_CORES_PER_SOCKET * NUM_THREADS_PER_CORE * socket_num)) * NUM_SOCKETS) / thread_mult + socket_num;

      cpu_set_t cpuset;
      CPU_ZERO(&cpuset);
      CPU_SET(cpu, &cpuset);

      /* int rc = */ pthread_setaffinity_np(thr, sizeof(cpu_set_t), &cpuset);
    }

    static
    size_t num_tries(size_t thisthr, size_t thatthr)
    {
      size_t thiscore = thisthr & thread_mask;
      size_t thatcore = thatthr & thread_mask;

      if (thiscore == thatcore) return 6;

      if ((thiscore < NUM_CORES_PER_SOCKET) == (thatcore < NUM_CORES_PER_SOCKET)) return 4;

      return 3;
    }

    static constexpr
    size_t max_concurrency()
    {
      return NUM_SOCKETS * NUM_CORES_PER_SOCKET * NUM_THREADS_PER_CORE;
    }
  };


  template <size_t NSOCKETS, size_t NCORES, size_t NTHREADS>
  struct power_arch
  {
    enum
    {
      NUM_SOCKETS          = NSOCKETS,
      NUM_CORES_PER_SOCKET = NCORES,
      NUM_THREADS_PER_CORE = NTHREADS,
      TOTAL_THREADS        = NUM_SOCKETS * NUM_CORES_PER_SOCKET * NUM_THREADS_PER_CORE
    };

    static
    void set_threadinfo_info(size_t numthreads)
    {
      assert(numthreads <= TOTAL_THREADS);

      thread_mult = std::min<size_t>(TOTAL_THREADS / numthreads, NUM_THREADS_PER_CORE);

      thread_mask = NUM_THREADS_PER_CORE / thread_mult;
      assert(thread_mask && (thread_mask & (thread_mask-1)) == 0);

      thread_mask = ~(thread_mask-1);
    }

    static
    void bind_to_core(pthread_t thr, size_t num)
    {
      num = num * NUM_THREADS_PER_CORE;
      num = num / TOTAL_THREADS + num % TOTAL_THREADS;

      cpu_set_t cpuset;
      CPU_ZERO(&cpuset);
      CPU_SET(num, &cpuset);

      /* int rc = */ pthread_setaffinity_np(thr, sizeof(cpu_set_t), &cpuset);
    }

    static
    size_t num_tries(size_t thisthr, size_t thatthr)
    {
      size_t thiscore = thisthr & thread_mask;
      size_t thatcore = thatthr & thread_mask;

      if (thiscore == thatcore) return 6;

      if ((thiscore < NUM_CORES_PER_SOCKET) == (thatcore < NUM_CORES_PER_SOCKET)) return 4;

      return 2;
    }

    static constexpr
    size_t max_concurrency()
    {
      return NUM_SOCKETS * NUM_CORES_PER_SOCKET * NUM_THREADS_PER_CORE;
    }
  };
#endif /* exclude non-linux platforms */
}
#endif /* _UCL_ARCHMODEL_H */
