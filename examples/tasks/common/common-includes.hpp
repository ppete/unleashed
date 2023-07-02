
#ifndef COMMON_INCLUDES_HPP
#define COMMON_INCLUDES_HPP 1

#include <sstream>

#if OMP_VERSION
#include <omp.h>

//~ #include "omp-continue.hpp"
#include "ucl/task-continue.hpp"
#endif /* OMP_VERSION */

#if WOMP_VERSION
// close to the original version in the Barcelona OpenMP Testing Suite (BOTS).
#include <omp.h>
#endif /* WOMP_VERSION */

#if TBB_VERSION
#include <tbb/task_scheduler_init.h>
#include <tbb/task_group.h>

#include "simple-reducer.hpp"  // most simple reducer, solely for benchmarking
#endif /* TBB_VERSION */

#if CILK_VERSION
#include <cstdio>
#include <cilk/cilk.h>
#include <cilk/reducer_opadd.h>
#include <cilk/cilk_api.h>
#endif /* CILK_VERSION */

#if QTHREADS_VERSION
#include <qthread/qthread.hpp>
#include <qthread/sinc.h>
#include "../common/qthreads.hpp"
#endif /* QTHREADS_VERSION */

#if QTHREADS_VERSION_ORIGINAL
#include <qthread/qthread.hpp>
#include <qthread/sinc.h>
#include "qthreads.hpp"
#endif /* QTHREADS_VERSION_ORIGINAL */

#if UCL_VERSION
  #define UCL_RUNTIME_DATA 0

  #include "ucl/task.hpp"
  #include "ucl/task-pool-x.hpp"
  #include "ucl/task-continue.hpp"
#endif /* UCL_VERSION */

#ifndef NUMTHREADS
#define NUMTHREADS 20
#endif /* NUMTHREADS */

namespace aux
{
  template <class T, class S>
  static inline
  T as(const S& t)
  {
    T                 res;
    std::stringstream str;

    str << t;
    str >> res;

    return res;
  }

  template <class T>
  static inline
  void unused(const T&) {}
}

#if CILK_VERSION

static inline
void set_cilk_workers(int n)
{
  assert(n <= 9999);

  char str[5];

  sprintf(str, "%d", n);

  bool failure = __cilkrts_set_param("nworkers", str) != 0;
  assert(!failure), aux::unused(failure);
}

static inline
void cilk_init(int workers, std::string stacksz)
{
  bool stackset = __cilkrts_set_param("stack size", stacksz.c_str()) == 0;
  assert(stackset), aux::unused(stackset);

  set_cilk_workers(workers);
}

#endif /* CILK_VERSION */


#endif /* COMMON_INCLUDES_HPP */
