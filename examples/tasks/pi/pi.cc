/**
 * Computes Pi through Numeric Integration
 *
 * The numeric integral code is based on Jim Lambers lecture notes on
 * Adaptive Quadrature.
 * http://www.math.usm.edu/lambers/mat460/fall09/lecture30
 *
 * Implementer: Peter Pirkelbauer (UAB)
 */

/**
 * This program is part of the Blaze-Task Test Suite
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
 * 3. Neither the name of the copyright holder nor the names of its contributors
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


#include <iostream>
#include <iomanip>
#include <chrono>
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <list>
#include <vector>

#define PRINT_STATS 0
#define WITH_HISTOGRAM 0

#ifndef NUMTHREADS
#define NUMTHREADS (20)
#endif /* NUMTHREADS */

#ifndef PROBLEM_SIZE
#define PROBLEM_SIZE (1e-17)
#endif /* PROBLEM_SIZE */


#if OMP_VERSION
#include <omp.h>
#endif /* OMP_VERSION */

#if TBB_VERSION
#include <mutex>
#include <tbb/task_scheduler_init.h>
#include <tbb/task_group.h>
#endif /* TBB_VERSION */

#if BLAZE_VERSION
  // NOTE: include archmodel.hpp and typedef arch_model to target system
  //       to make number of work-stealing attempts sensitive to
  //       thief-victim cache hierachy. Mileage varies depending on
  //       benchmark.
  //~ #include "archmodel.hpp"

  //~ typedef uab::power_arch<2, 20, 4>   arch_model; // power9 dual socket
  //~ typedef uab::power_arch<2, 10, 8>   arch_model; // power8 dual socket
  //~ typedef uab::intel_arch<2, 10, 2> arch_model;    // intel dual socket

  #include "tasks.hpp"
#endif /* BLAZE_VERSION */

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

#include "atomicutil.hpp"


typedef long double pi_type;

static const pi_type EPSILON    = PROBLEM_SIZE; ///< controls how many tasks will be needed

#if WITH_HISTOGRAM
template <class D>
struct histogram
{
  D                   bucketsz;
  std::vector<size_t> hstgrm;

  histogram(D lo, D hi, size_t buckets)
  : bucketsz((hi-lo)/buckets), hstgrm(buckets+1, 0)
  {
    assert(hstgrm.size() > buckets);
  }

  void add(D pos) { ++hstgrm.at(pos/bucketsz); }

  void add(D lo, D step)
  {
    add(lo);
    add(lo+step);
  }

  void print()
  {
    std::cout << hstgrm.at(0);

    for (size_t i = 1; i < hstgrm.size(); ++i)
    {
      std::cout << " . " << hstgrm[i];
    }

    std::cout << std::endl;
  }
};

histogram<long double> hist(0.0, 1.0, 100);

#endif /* WITH_HISTOGRAM */

/// number of segments used to compute the result
static size_t        splits = 1;

struct pi_formula
{
  template <class D>
  D operator()(D x)
  {
    return sqrt(1-x*x);
  }
};

template <class F, class D>
auto rectangular(F f, D low, D step, size_t iter) -> D
{
  D x = low + (iter+D(.5))*step;
  return f(x) * step;
}

template <class F, class D>
auto trapezoidal(F f, D low, D step, size_t iter) -> D
{
  D xa = low + iter*step;

  return (step/2) * (f(xa) + f(xa+step));
}

template <class F, class D>
auto method(F f, D low, D step, size_t iter) -> D
{
  // return rectangular(f, low, step, iter);
  return trapezoidal(f, low, step, iter);
}


#if NOTASK

static const size_t  STEPS      = 1982;

template <class D, class F>
auto integrate(F f, D low, D hi, D, size_t steps = STEPS) -> D
{
  D res  = 0.0;
  D step = (hi-low) / steps;

  #pragma omp parallel for reduction (+:res)
  for (size_t i = 0; i < steps; ++i)
  {
    res += method(f, low, step, i);
  }

  return res;
}

#endif /* NOTASK */

template <class D>
struct integration_task
{
  D low;
  D step;
  D res;
};


#if OMP_VERSION

static inline
size_t thread_num()
{
  return omp_get_thread_num();
}

template <class F, class D>
auto integrate_adaptive( F f, D eps, integration_task<D> task,
                         uab::aligned_type<std::pair<D, size_t>, CACHELINESZ>* s
                       ) -> void // std::pair<D, size_t>
{
  typedef integration_task<D> integration_task;

  size_t thrnum = thread_num();
  D      dif;
  D      tol;
  D      a;

  do
  {
    D    halfstep = task.step / 2;
    D    a1       = method(f, task.low,          halfstep, 0);
    D    a2       = method(f, task.low+halfstep, halfstep, 0);

    a        = a1+a2;
    dif = task.res > a ? task.res - a : a - task.res;
    tol = 3 * task.step * eps;

    if (dif >= tol)
    {
      #pragma omp task
      integrate_adaptive(f, eps, integration_task{task.low,          halfstep, a1}, s);

      ++s[thrnum].val.second;

      // run on same thread
      task = integration_task{task.low+halfstep, halfstep, a2};
    }
  } while (dif >= tol);

  s[thrnum].val.first += a;

  //~ hist.add(task.low, task.step);
}


template <class D, class F>
auto integrate_adaptive(F f, D lo, D hi, D eps) -> D
{
  omp_set_num_threads(NUMTHREADS);

  uab::aligned_type<std::pair<D, size_t>, CACHELINESZ> s[NUMTHREADS];
  D                                                    step = hi-lo;
  D                                                    res  = D();
  size_t                                               ctr  = 1;

  #pragma omp parallel firstprivate(f, lo, hi, eps, step) shared(res,ctr,s)
  {
    size_t thrnum = thread_num();

    s[thrnum].val = std::pair<D, size_t>(D(), 0);

    #pragma omp single
    {
      #pragma omp taskgroup
      {
        integrate_adaptive(f, eps, integration_task<D>{lo, step, method(f, lo, step, 0)}, s);
      }
    }

    #pragma omp atomic
    res += s[thrnum].val.first;

    #pragma omp atomic
    ctr += s[thrnum].val.second;
  }

  splits = ctr;

  //~ hist.print();

#if PRINT_STATS
  for (size_t i = 0; i < NUMTHREADS; ++i)
    std::cout << i << ": " << s[i].val.second << std::endl;
#endif

  return res;
}
#endif /* OMP_VERSION */

#if TBB_VERSION

template <class D> struct globals
{
  static uab::aligned_atomic_type<D, CACHELINESZ>      result;
  static uab::aligned_atomic_type<size_t, CACHELINESZ> splits;

  static void add(D d)
  {
    D v = result.val.load(std::memory_order_relaxed);

    while (!result.val.compare_exchange_strong(v, v+d, std::memory_order_relaxed, std::memory_order_relaxed));
  }
};

template <class D>
uab::aligned_atomic_type<D, CACHELINESZ> globals<D>::result(0);


template <class D>
uab::aligned_atomic_type<size_t, CACHELINESZ> globals<D>::splits(1);


template <class G, class F, class D>
auto integrate_adaptive(G& taskgroup, F f, D eps, integration_task<D> task) -> void // std::pair<D, size_t>
{
  typedef integration_task<D> integration_task;

  D dif;
  D tol;
  D a;

  do
  {
    D halfstep = task.step / 2;
    D a1       = method(f, task.low,          halfstep, 0);
    D a2       = method(f, task.low+halfstep, halfstep, 0);

    a   = a1+a2;
    dif = task.res > a ? task.res - a : a - task.res;
    tol = 3*task.step*eps;

    if (dif >= tol)
    {
      taskgroup.run([&taskgroup,f,eps,task,halfstep,a1]()->void
                     { integrate_adaptive(taskgroup, f, eps, integration_task{task.low, halfstep, a1});
                     }
                   );

#ifndef ZEROGLOBALS
      globals<D>::splits.val.fetch_add(1, std::memory_order_relaxed);
#endif /* !ZEROGLOBALS */

      // run on same thread
      task = integration_task{task.low+halfstep, halfstep, a2};
    }
  } while (dif >= tol);

#ifndef ZEROGLOBALS
      globals<D>::add(a);
#endif /* !ZEROGLOBALS */
}


template <class D, class F>
auto integrate_adaptive(F f, D lo, D hi, D eps) -> D
{
  tbb::task_scheduler_init init(NUMTHREADS);
  tbb::task_group          g;
  D                        step = hi-lo;

  integrate_adaptive(g, f, eps, integration_task<D>{lo, step, method(f, lo, step, 0)});

  g.wait();
  splits = globals<D>::splits.val.load(std::memory_order_relaxed);
  return globals<D>::result.val.load(std::memory_order_relaxed);
}

#endif /* TBB_VERSION */

#if CILK_VERSION

void set_cilk_workers(int n)
{
  assert(n <= 9999);

  char str[5];

  sprintf(str, "%d", n);

  bool success = __cilkrts_set_param("nworkers", str) != 0;
  assert(success);
}

template <class F, class D>
auto integrate_adaptive(F f, D eps, integration_task<D> task, cilk::reducer_opadd<D>& sum) -> void
{
  typedef integration_task<D> integration_task;

  D dif;
  D tol;
  D res;

  do
  {
    D halfstep = task.step / 2;
    D a1       = method(f, task.low,          halfstep, 0);
    D a2       = method(f, task.low+halfstep, halfstep, 0);

    res = a1+a2;
    dif = task.res > res ? task.res - res : res - task.res;
    tol = 3*task.step*eps;

    if (dif >= tol)
    {
      cilk_spawn integrate_adaptive(f, eps, integration_task{task.low, halfstep, a1}, sum);

      // run on same thread
      task = integration_task{task.low+halfstep, halfstep, a2};
    }
  } while (dif >= tol);

  sum += res;
}

template <class D, class F>
auto integrate_adaptive(F f, D lo, D hi, D eps) -> D
{
  set_cilk_workers(NUMTHREADS);

  D                      step = hi-lo;
  cilk::reducer_opadd<D> sum;

  integrate_adaptive(f, eps, integration_task<D>{lo, step, method(f, lo, step, 0)}, sum);

  return sum.get_value();
}

#endif /* CILK_VERSION */

#if QTHREADS_VERSION

template <class F, class D>
struct qtask
{
  F f;
  D eps;
  D low;
  D step;
  D res;
  qt_sinc_t*       sinc;
};

template <class F, class D>
aligned_t qthreads_integrate(void* qtsk)
{
  typedef qtask<F,D> qtask;

  qtask& task = *reinterpret_cast<qtask*>(qtsk);

  D dif;
  D tol;
  D res;

  do
  {
    D halfstep = task.step / 2;
    D a1       = method(task.f, task.low,          halfstep, 0);
    D a2       = method(task.f, task.low+halfstep, halfstep, 0);

    res = a1+a2;
    dif = task.res > res ? task.res - res : res - task.res;
    tol = 3*task.step*task.eps;

    if (dif >= tol)
    {
      qt_sinc_expect(task.sinc, 1);
      qthread_fork( qthreads_integrate<F,D>,
                    new qtask{ task.f, task.eps, task.low, halfstep, a1, task.sinc},
                    nullptr
                  );

      // run on same thread
      task = qtask{ task.f, task.eps, task.low+halfstep, halfstep, a2, task.sinc};
    }
  } while (dif >= tol);

  qt_sinc_submit(task.sinc, &res);
  delete &task;
  return 0;
}

template <class D>
void reduce(void* target, const void* source)
{
  D*       tgt = reinterpret_cast<D*>(target);
  const D* src = reinterpret_cast<const D*>(source);

  *tgt += *src;
}

template <class D, class F>
auto integrate_adaptive(F f, D lo, D hi, D eps) -> D
{
  D          result = 0;
  D          step   = hi-lo;
  qt_sinc_t* sinc   = qt_sinc_create(sizeof(D), &result, reduce<D>, 1);

  qthreads_integrate<F,D>(new qtask<F,D>{f, eps, lo, step, method(f, lo, step, 0), sinc});

  qt_sinc_wait(sinc, &result);
  return result;
}

#endif /* QTHREADS_VERSION */



#if BLAZE_VERSION

template <class D, class F, class T>
struct adaptive_integral
{
  F        fun;
  D        eps;

  explicit
  adaptive_integral(F f, D e)
  : fun(f), eps(e)
  {}

  template <class pool_t>
  D operator()(pool_t& tasks, T task)
  {
    D   dif;
    D   tol;
    D   res;

    do
    {
      D halfstep   = task.step / 2;
      D a1         = method(fun, task.low,          halfstep, 0);
      D a2         = method(fun, task.low+halfstep, halfstep, 0);

      res = a1 + a2;
      dif = task.res > res ? task.res - res : res - task.res;
      tol = 3*task.step*eps;

      if (dif >= tol)
      {
        tasks.enq(T{task.low,          halfstep, a1});
        task =    T{task.low+halfstep, halfstep, a2};
      }
    } while (dif >= tol);

    return res;
  }
};

template <class F, class D>
D integrate_adaptive(F f, D lo, D hi, D eps)
{
  adaptive_integral<D, F, integration_task<D> > fun(f, eps);

  return uab::execute_tasks(NUMTHREADS, fun, integration_task<D>{ lo, hi-lo, eps });
}

#endif /* BLAZE_VERSION */

int main()
{
  typedef std::chrono::time_point<std::chrono::system_clock> time_point;

#if QTHREADS_VERSION
  init_qthreads(NUMTHREADS);
#endif /* QTHREADS_VERSION */

  time_point     starttime = std::chrono::system_clock::now();

  // executes loop in parallel
  //   and uses a reduction algorithm to combine all pi values that were
  //   computed across threads.
  pi_type        pi  = 4 * integrate_adaptive(pi_formula(), pi_type(0), pi_type(1), EPSILON);

  time_point     endtime = std::chrono::system_clock::now();
  int            elapsedtime = std::chrono::duration_cast<std::chrono::milliseconds>(endtime-starttime).count();

  std::cout << "pi   = " << std::setprecision(32) << pi << std::endl;
  std::cout << "time = " << elapsedtime << "ms" << std::endl;
  std::cout << splits << std::endl;
  std::cerr << elapsedtime << std::endl;
  return 0;
}
