
#ifndef _ARCHMODEL_H
#define _ARCHMODEL_H 1

#include <cassert>

//
// provides implementation of hardware models

namespace ucl
{
  size_t thread_mult;  // == NUM_THREADS_PER_CORE if hyperthreading is not needed
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
  };

#if !defined __CYGWIN__ && !defined OS_WIN32 && !defined __OpenBSD__ && !defined __sun
  template <size_t NSOCKETS, size_t NCORES, size_t NTHREADS>
  struct intel_arch
  {
    enum 
    { 
      NUM_SOCKETS          = NSOCKETS,
      NUM_CORES_PER_SOCKET = NCORES,
      NUM_THREADS_PER_CORE = NTHREADS,
      NUM_CPUS             = NUM_CORES_PER_SOCKET * NTHREADS,
      NUM_CORES            = NUM_CPUS * NTHREADS
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

      cpu = cpu * thread_mult;

      const size_t socket_num = size_t(cpu >= (NUM_CORES_PER_SOCKET * NUM_THREADS_PER_CORE));

      cpu = (cpu - (NUM_CORES_PER_SOCKET * NUM_THREADS_PER_CORE * socket_num)) * NUM_SOCKETS / thread_mult + socket_num;

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
  };


  template <size_t NSOCKETS, size_t NCORES, size_t NTHREADS>
  struct power_arch
  {
    enum 
    { 
      NUM_SOCKETS          = NSOCKETS,
      NUM_CORES_PER_SOCKET = NCORES,
      NUM_THREADS_PER_CORE = NTHREADS,
      NUM_CPUS             = NUM_CORES_PER_SOCKET * NTHREADS,
      NUM_CORES            = NUM_CPUS * NTHREADS
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
    void bind_to_core(pthread_t thr, size_t num)
    {
      static_assert((NTHREADS & (NTHREADS-1)) == 0, "assumed power of 2");

      //  0 -> 0,  1 ->  8, 2 -> 16, .., 19 -> 152
      // 20 -> 1, 21 ->  9,
      // 40 -> 2, 41 -> 10,

      num = num * thread_mult;

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
  };
#endif /* exclude non-linux platforms */
}
#endif /* _ARCHMODEL_H */
