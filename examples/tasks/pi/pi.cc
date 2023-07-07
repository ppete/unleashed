/**
 * Computes Pi through Numeric Integration
 *
 * The numeric integral code is based on Jim Lambers lecture notes on
 * Adaptive Quadrature.
 * http://www.math.usm.edu/lambers/mat460/fall09/lecture30
 *
 * Implementer: Peter Pirkelbauer
 */

/**
 * The Unleashed Concurrency Library's Task testing framework
 *
 * Copyright (c) 2019, LLNL
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

#include "../common/common-includes.hpp"

#include "ucl/atomicutil.hpp"

#define PRINT_STATS 0
#define WITH_HISTOGRAM 0

#ifndef PROBLEM_SIZE
#define PROBLEM_SIZE (1e-17)
#endif /* PROBLEM_SIZE */


using pi_type = long double;

struct counting_task_t
{
  counting_task_t(pi_type v = 0.0)
  : val(v), seg(1)
  {}

  counting_task_t& operator+=(const counting_task_t& that)
  {
    this->val += that.val;
    this->seg += that.seg;
    return *this;
  }

  friend
  counting_task_t operator+(const counting_task_t& lhs, const counting_task_t& rhs)
  {
    counting_task_t tmp(lhs);

    return tmp+=rhs;
  }

  pi_type val;
  size_t  seg;
};

typedef pi_type                    task_result_type;
//~ typedef counting_task_t            task_result_type;


namespace
{
inline
pi_type value(counting_task_t t)   { return t.val; }

inline
size_t segments(counting_task_t t) { return t.seg; }

inline
pi_type value(pi_type t)   { return t; }

inline
size_t segments(pi_type)   { return 0; }
}



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
struct pi_formula
{
  template <class D>
  D operator()(D x)
  {
    return sqrt(1-x*x);
  }

#if __PGI
  pi_formula(const pi_formula&) {}

  pi_formula() {}
#endif
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

template <class D>
struct integration_task
{
  D low;
  D step;
  D res;

#if UCL_VERSION
  integration_task(const integration_task&) = delete;
  integration_task& operator=(const integration_task&) = delete;

  integration_task(D l, D s, D r)
  : low(l), step(s), res(r)
  {}

  integration_task()
  : integration_task(D(), D(), D())
  {}

  integration_task(integration_task&& other)
  : low(other.low), step(other.step), res(other.res)
  {}

  integration_task& operator=(integration_task&& other)
  {
    low = other.low;
    step = other.step;
    res = other.res;
    return *this;
  }
#elif __PGI
  integration_task(const integration_task& other)
  : low(other.low), step(other.step), res(other.res)
  {}

  integration_task(D l, D s, D r)
  : low(l), step(s), res(r)
  {}
#endif /* UCL_VERSION */
};


#if OMP_VERSION

pi_type partialresult;
#pragma omp threadprivate(partialresult)

template <class F, class D>
auto integrate_adaptive(F f, D eps, integration_task<D> task) -> void // std::pair<D, size_t>
{
  typedef integration_task<D> integration_task;

  D      dif;
  D      tol;
  D      a;

  do
  {
    D    halfstep = task.step / 2;
    D    a1       = method(f, task.low,          halfstep, 0);
    D    a2       = method(f, task.low+halfstep, halfstep, 0);

    a   = a1+a2;
    dif = task.res > a ? task.res - a : a - task.res;
    tol = 3 * task.step * eps;

    if (dif >= tol)
    {
      #pragma omp task firstprivate(f, eps, task, halfstep, a1)
      integrate_adaptive(f, eps, integration_task{task.low, halfstep, a1});

      // run on same thread
      task = integration_task{task.low+halfstep, halfstep, a2};
    }
  } while (dif >= tol);

  partialresult += a;

  //~ hist.add(task.low, task.step);
}


template <class D, class F>
auto integrate_adaptive(F f, D lo, D hi, size_t numthreads, D eps) -> task_result_type
{
  D step = hi-lo;
  D res  = D();

  #pragma omp parallel num_threads(numthreads) firstprivate(f, lo, eps, step) shared(res)
  {
    partialresult = D();

    #pragma omp single
    #pragma omp taskgroup
    {
      integrate_adaptive(f, eps, integration_task<D>{lo, step, method(f, lo, step, 0)});
    }

    #pragma omp atomic
    res += partialresult;
  }

  //~ hist.print();
  return res;
}
#endif /* OMP_VERSION */

#if WOMP_VERSION || SEQ_VERSION

pi_type partialresult;
#pragma omp threadprivate(partialresult)

template <class F, class D>
auto integrate_adaptive(F f, D eps, integration_task<D> task) -> void // std::pair<D, size_t>
{
  typedef integration_task<D> integration_task;

  D      dif;
  D      tol;
  D      a;

  do
  {
    D    halfstep = task.step / 2;
    D    a1       = method(f, task.low,          halfstep, 0);
    D    a2       = method(f, task.low+halfstep, halfstep, 0);

    a   = a1+a2;
    dif = task.res > a ? task.res - a : a - task.res;
    tol = 3 * task.step * eps;

    if (dif >= tol)
    {
      #pragma omp task
      integrate_adaptive(f, eps, integration_task{task.low, halfstep, a1});

      // run on same thread
      task = integration_task{task.low+halfstep, halfstep, a2};
    }
  } while (dif >= tol);

  #pragma omp taskwait

  partialresult += a;

  //~ hist.add(task.low, task.step);
}


template <class D, class F>
auto integrate_adaptive(F f, D lo, D hi, size_t numthreads, D eps) -> task_result_type
{
  D step = hi-lo;
  D res  = D();

  #pragma omp parallel num_threads(numthreads) firstprivate(f, lo, eps, step) shared(res)
  {
    partialresult = D();

    #pragma omp single
    #pragma omp taskgroup
    {
      integrate_adaptive(f, eps, integration_task<D>{lo, step, method(f, lo, step, 0)});
    }

    #pragma omp atomic
    res += partialresult;
  }

  //~ hist.print();
  return res;
}
#endif /* OMP_VERSION */


#if TBB_VERSION

template <class G, class F, class D>
void integrate_adaptive( G& taskgroup,
                         F f,
                         D eps,
                         integration_task<D> task,
                         ucl::simple_reducer<D>& reducer
                       )
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
      taskgroup.run( [&taskgroup,f,eps,task,halfstep,a1,&reducer]()->void
                     {
                       integrate_adaptive( taskgroup,
                                           f,
                                           eps,
                                           integration_task{task.low, halfstep, a1},
                                           reducer
                                         );
                     }
                   );

      // run on same thread
      task = integration_task{task.low+halfstep, halfstep, a2};
    }
  } while (dif >= tol);

  reducer += a;
}


template <class D, class F>
auto integrate_adaptive(F f, D lo, D hi, size_t numthreads, D eps) -> task_result_type
{
  TBB_INIT(numthreads);
  tbb::task_group                       g;
  ucl::simple_reducer<task_result_type> reducer;
  D                                     step = hi-lo;

  integrate_adaptive(g, f, eps, integration_task<D>{lo, step, method(f, lo, step, 0)}, reducer);

  g.wait();
  return reducer.get_value();
}

#endif /* TBB_VERSION */

#if CILK_VERSION

void sum_init(void* sum) { *static_cast<pi_type*>(sum) = 0.0; }
void sum_plus(void* lhs, void* rhs) { *static_cast<pi_type*>(lhs) += *static_cast<pi_type*>(rhs); }

pi_type cilk_reducer(sum_init, sum_plus) sum(0);

void sum_add(pi_type addend)
{
  sum += addend;
}

pi_type sum_result()
{
  return sum;
}

template <class F, class D>
auto integrate_adaptive( F f, D eps, integration_task<D> task
                       ) -> void
{
  using integration_task = integration_task<D>;

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
      cilk_spawn integrate_adaptive(f, eps, integration_task{task.low, halfstep, a1});

      // run on same thread
      task = integration_task{task.low+halfstep, halfstep, a2};
    }
  } while (dif >= tol);

  sum_add(res);
}

template <class D, class F>
auto integrate_adaptive(F f, D lo, D hi, size_t /*numthreads*/, D eps) -> task_result_type
{
  D       step = hi-lo;

  integrate_adaptive(f, eps, integration_task<D>{lo, step, method(f, lo, step, 0)});
  return sum_result();
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
auto integrate_adaptive(F f, D lo, D hi, size_t /*numthreads*/, D eps) -> D
{
  // numthreads is set at start up

  D          result = 0;
  D          step   = hi-lo;
  qt_sinc_t* sinc   = qt_sinc_create(sizeof(D), &result, reduce<D>, 1);

  qthreads_integrate<F,D>(new qtask<F,D>{f, eps, lo, step, method(f, lo, step, 0), sinc});

  qt_sinc_wait(sinc, &result);
  return result;
}

#endif /* QTHREADS_VERSION */



#if UCL_VERSION

template <class D, class F, class T>
struct adaptive_integral
{
  //~ typedef D result_type;
  // typedef std::pair<D, size_t> result_type;

  F        fun;
  D        eps;

  explicit
  adaptive_integral(F f, D e)
  : fun(f), eps(e)
  {}

  template <class pool_t>
  task_result_type operator()(pool_t& tasks, T task)
  {
    D dif;
    D tol;
    D res;

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
task_result_type integrate_adaptive(F f, D lo, D hi, size_t numthreads, D eps)
{
  adaptive_integral<D, F, integration_task<D> > fun(f, eps);

  return ucl::execute_tasks_x(numthreads, fun, integration_task<D>{ lo, hi-lo, eps });
}

#endif /* UCL_VERSION */

int main(int argc, char** args)
{
  typedef std::chrono::time_point<std::chrono::system_clock> time_point;

  size_t  num_threads = NUMTHREADS;
  pi_type tolerance   = EPSILON;

  if (argc > 1) num_threads = aux::as<size_t>(*(args+1));
  if (argc > 2) tolerance   = aux::as<pi_type>(*(args+2));

  std::cout << num_threads << " << threads - tolerance >> " << tolerance
            << std::endl;

#if QTHREADS_VERSION
  init_qthreads(num_threads);
#endif /* QTHREADS_VERSION */

  time_point     starttime = std::chrono::system_clock::now();

  // executes loop in parallel
  //   and uses a reduction algorithm to combine all pi values that were
  //   computed across threads.
  task_result_type pi  = integrate_adaptive( pi_formula{},
                                             pi_type(0),
                                             pi_type(1),
                                             num_threads,
                                             tolerance
                                           );

  time_point     endtime = std::chrono::system_clock::now();
  int            elapsedtime = std::chrono::duration_cast<std::chrono::milliseconds>(endtime-starttime).count();

  std::cout << "pi   = " << std::setprecision(32) << (4*value(pi)) << std::endl;
  std::cout << "segs = " << segments(pi) << std::endl;
  std::cout << "time = " << elapsedtime << "ms" << std::endl;

  std::cerr << elapsedtime << std::endl;
  return 0;
}
