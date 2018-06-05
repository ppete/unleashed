// openmp:     g++ -std=c++11 -Wall -Wextra -pedantic -O2 -fopenmp -DOMP_VERSION=1 pi.cc -o /tmp/pi.bin
// blaze:      g++ -std=c++11 -Wall -Wextra -pedantic -O2 -pthread -DBLAZE_VERSION=1 pi.cc -o /tmp/pi.bin
// tbb:        g++ -std=c++11 -Wall -Wextra -pedantic -O2 -pthread -DTBB_VERSION=1 -I$TBB_HOME/include -L$TBB_HOME/lib -ltbb pi.cc -o /tmp/pi.bin
// sequential? g++ -std=c++11 -Wall -Wextra -pedantic -O2 -pthread pi.cc -o /tmp/pi.bin

#include <iostream>
#include <iomanip>
#include <chrono>
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <thread>
#include <list>
#include <vector>

#define PRINT_STATS 1

#ifndef NUMTHREADS
#define NUMTHREADS (20)
#endif /* NUMTHREADS */


#if OMP_VERSION
#include <omp.h>
#endif

#if TBB_VERSION
#include <mutex>
#include <tbb/task_scheduler_init.h>
#include <tbb/task_group.h>
#endif

#if BLAZE_VERSION

#include "tasks.hpp"

#endif

#include "atomicutil.hpp"


#if 1 /* WITH_HISTOGRAM */
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
#endif /* WITH_HISTOGRAM */

typedef long double pi_type;

static const pi_type EPSILON    = 1e-17; ///< controls how many tasks will be needed

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

//~ histogram<long double> hist(0.0, 1.0, 100);

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
  D               step = hi-lo;
  tbb::task_group g;

  integrate_adaptive(g, f, eps, integration_task<D>{lo, step, method(f, lo, step, 0)});

  g.wait();
  splits = globals<D>::splits.val.load(std::memory_order_relaxed);
  return globals<D>::result.val.load(std::memory_order_relaxed);
}

#endif /* TBB_VERSION */


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

  time_point     starttime = std::chrono::system_clock::now();

  // executes loop in parallel
  //   and uses a reduction algorithm to combine all pi values that were
  //   computed across threads.
  #if OMP_VERSION
  omp_set_num_threads(NUMTHREADS);

  pi_type        pi  = 4 * integrate_adaptive(pi_formula(), pi_type(0), pi_type(1), EPSILON);
  #endif

  #if TBB_VERSION
  tbb::task_scheduler_init init(NUMTHREADS);

  pi_type        pi  = 4 * integrate_adaptive(pi_formula(), pi_type(0), pi_type(1), EPSILON);
  #endif

  #if BLAZE_VERSION
  pi_type        pi  = 4 * integrate_adaptive(pi_formula(), pi_type(0), pi_type(1), EPSILON);
  #endif

  time_point     endtime = std::chrono::system_clock::now();
  int            elapsedtime = std::chrono::duration_cast<std::chrono::milliseconds>(endtime-starttime).count();

  std::cout << "pi   = " << std::setprecision(32) << pi << std::endl;
  std::cout << "time = " << elapsedtime << "ms" << std::endl;
  std::cout << splits << std::endl;
  std::cerr << elapsedtime << std::endl;
  return 0;
}
