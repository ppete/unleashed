
#ifndef COMMON_INCLUDES_HPP
#define COMMON_INCLUDES_HPP 1

#ifndef NUMTHREADS
#define NUMTHREADS (20)
#endif /* NUMTHREADS */

#if OMP_VERSION
#include <omp.h>
#endif /* OMP_VERSION */

#if BOTS_VERSION
// close to the original version in the Barcelona OpenMP Testing Suite (BOTS).
#include <omp.h>
#endif /* BOTS_VERSION */

#if TBB_VERSION
#include <tbb/task_scheduler_init.h>
#include <tbb/task_group.h>

#include "../common/simple-reducer.hpp"  // most simple reducer, solely for benchmarking
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
#include "../common/qthreads.hpp"
#endif /* QTHREADS_VERSION_ORIGINAL */

#if BLAZE_VERSION
  // NOTE: include archmodel.hpp and typedef arch_model to target system
  //       to make number of work-stealing attempts sensitive to
  //       thief-victim cache hierarchy. Mileage varies depending on
  //       benchmark.
  // #include "archmodel.hpp"

  // typedef uab::power_arch<2, 20, 4>   arch_model; // power9 dual socket
  // typedef uab::power_arch<2, 10, 8>   arch_model; // power8 dual socket
  // typedef uab::intel_arch<2, 10, 2> arch_model;   // intel dual socket

  #include "tasks.hpp"
#endif /* BLAZE_VERSION */

#endif /* COMMON_INCLUDES_HPP */
