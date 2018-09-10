/**
 * A simple task based implementation to compute fibonacci numbers
 *
 * Implementer: Peter Pirkelbauer (UAB) - 2018
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
#include <thread>
#include <list>
#include <vector>

#ifndef NUMTHREADS
#define NUMTHREADS (20)
#endif /* NUMTHREADS */

#ifndef PROBLEM_SIZE
#define PROBLEM_SIZE (40)
#endif /* PROBLEM_SIZE */

#define PRINT_STATS 0

#if OMP_VERSION
#include <omp.h>
#endif

#if TBB_VERSION
#include <tbb/task_scheduler_init.h>
#include <tbb/task_group.h>
#endif

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
#include <sstream>
#include <qthread/qthread.hpp>
#include <qthread/sinc.h>
#endif /* QTHREADS_VERSION */

#include "atomicutil.hpp"


typedef long long fib_type;

#if 0 /* WITH_HISTOGRAM */
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


template <class I>
struct fibonacci_task
{
  I num;
};


#if OMP_VERSION

//~ histogram<long double> hist(0.0, 1.0, 100);

static inline
size_t thread_num()
{
  return omp_get_thread_num();
}

template <class I>
void fib_task(fibonacci_task<I> task, uab::aligned_type<std::pair<I, size_t>, CACHELINESZ>* s)
{
  typedef fibonacci_task<I> fibonacci_task;

  size_t thrnum = thread_num();

  while (task.num > 1)
  {
    #pragma omp task
    fib_task(fibonacci_task{task.num-2}, s);

    task = fibonacci_task{task.num-1};
  }

  s[thrnum].val.first  += task.num;
  s[thrnum].val.second += 1;
}

template <class I>
auto fib_task(I num) -> std::pair<I, size_t>
{
  uab::aligned_type<std::pair<I, size_t>, CACHELINESZ> s[NUMTHREADS];
  I                                                    res  = 0;
  size_t                                               ctr  = 1;

  omp_set_num_threads(NUMTHREADS);

  #pragma omp parallel firstprivate(num) shared(res,ctr,s)
  {
    size_t thrnum = thread_num();

    assert(thrnum < NUMTHREADS);
    s[thrnum] = std::make_pair(I(0), size_t(0));

    #pragma omp single
    {
      #pragma omp taskgroup
      {
        fib_task(fibonacci_task<I>{num}, s);
      }
    }

    #pragma omp atomic
    res += s[thrnum].val.first;

    #pragma omp atomic
    ctr += s[thrnum].val.second;
  }

  //~ hist.print();

#if PRINT_STATS
  for (size_t i = 0; i < NUMTHREADS; ++i)
    std::cout << i << ": " << s[i].val.second << std::endl;
#endif

  return std::make_pair(res, ctr);
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


template <class G, class I>
auto compute_fib(G& taskgroup, I task) -> void // std::pair<D, size_t>
{
  while (task > 1)
  {
    taskgroup.run( [&taskgroup, task]()->void
                   { compute_fib(taskgroup, task-2);
                   }
                 );
    task =    task-1;
  }
  // return task;
}

template <class I>
auto fib_task(I num) -> std::pair<I,size_t>
{
  tbb::task_scheduler_init init(NUMTHREADS);
  tbb::task_group          g;

  compute_fib(g, num);

  g.wait();
  return std::pair<I,size_t>();
}

#endif /* TBB_VERSION */


#if BLAZE_VERSION

template <class I>
struct compute_fib
{
  template <class T, class A>
  I operator()(uab::pool<T,A>& tasks, T task)
  {
    typedef fibonacci_task<I> fibonacci_task;

    while (task.num > 1)
    {
      tasks.enq(fibonacci_task{task.num-2});
      task =    fibonacci_task{task.num-1};
    }

    return task.num;
  }
};

template <class I>
std::pair<I, size_t> fib_task(I num)
{
  typedef fibonacci_task<I> fibonacci_task;

  compute_fib<I> fun;

  return std::make_pair(uab::execute_tasks(NUMTHREADS, fun, fibonacci_task{ num }), 0);
}

#endif /* BLAZE_VERSION */

#if CILK_VERSION

void set_cilk_workers(int n)
{
  assert(n <= 9999);

  char str[5];

  sprintf(str, "%d", n);

  bool success = __cilkrts_set_param("nworkers", str) != 0;
  assert(success);
}


template <class I>
void compute_fib(fibonacci_task<I> task, cilk::reducer_opadd<I>& sum)
{
  typedef fibonacci_task<I> fibonacci_task;

  while (task.num > 1)
  {
    cilk_spawn compute_fib(fibonacci_task{task.num-2}, sum);

    task =    fibonacci_task{task.num-1};
  }

  sum += task.num;
}

template <class I>
std::pair<I, size_t> fib_task(I num)
{
  typedef fibonacci_task<I> fibonacci_task;

  set_cilk_workers(NUMTHREADS);

  cilk::reducer_opadd<I> sum;

  compute_fib<I>(fibonacci_task{ num }, sum);
  return std::make_pair(sum.get_value(), 0);
}

#endif /* CILK_VERSION */


#if QTHREADS_VERSION

template <class I>
struct qtask
{
  I          num;
  qt_sinc_t* sinc;
};

void init_qthreads()
{
  std::stringstream str;

  str << "QTHREAD_HWPAR=" << NUMTHREADS;

  char* envset = new char[str.str().size()+1];

  memcpy(envset, str.str().c_str(), str.str().size()+1);
  putenv(envset);

  qthread_initialize();
}


template <class I>
aligned_t compute_fib(void* qtsk)
{
  qtask<I>& task = *reinterpret_cast<qtask<I>*>(qtsk);

  while (task.num > 1)
  {
    qtask<I> tsk2 = {task.num-2, task.sinc};

    qt_sinc_expect(task.sinc, 1);
    qthread_fork_copyargs( compute_fib<I>, &tsk2, sizeof(qtask<I>), nullptr );

    task.num = task.num-1;
  }

  qt_sinc_submit(task.sinc, &task.num);
  return 0;
}

template <class I>
void reduce(void* target, const void* source)
{
  I*       tgt = reinterpret_cast<I*>(target);
  const I* src = reinterpret_cast<const I*>(source);

  *tgt += *src;
}


template <class I>
std::pair<I, size_t> fib_task(I num)
{
  I          result;
  qt_sinc_t* sinc   = qt_sinc_create(sizeof(I), &result, reduce<I>, 1);
  qtask<I>   task   = { num, sinc };

  qthread_fork_copyargs(compute_fib<I>, &task, sizeof(qtask<I>), nullptr);
  qt_sinc_wait(sinc, &result);
  return std::make_pair(result, 0);
}

#endif /* QTHREADS_VERSION */



int main()
{
  typedef std::chrono::time_point<std::chrono::system_clock> time_point;

#if QTHREADS_VERSION
  init_qthreads();
#endif

  time_point     starttime = std::chrono::system_clock::now();

  // executes loop in parallel
  //   and uses a reduction algorithm to combine all pi values that were
  //   computed across threads.
  std::pair<fib_type, size_t> num  = fib_task(PROBLEM_SIZE);

  time_point     endtime = std::chrono::system_clock::now();
  int            elapsedtime = std::chrono::duration_cast<std::chrono::milliseconds>(endtime-starttime).count();

  std::cout << "num    = " << PROBLEM_SIZE << std::endl;
  std::cout << "fib    = " << num.first    << std::endl;
  std::cout << "leaves = " << num.second   << std::endl;
  std::cout << "time   = " << elapsedtime  << "ms" << std::endl;
  std::cerr << elapsedtime << std::endl;
  return 0;
}
