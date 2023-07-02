/**
 * A simple task based implementation to compute fibonacci numbers
 *
 * Implementer: Peter Pirkelbauer (LLNL)
 *
 * Implementer: The code under QTHREADS_ORIGINAL was taken from the qthreads
 *    distribution. Qthreads are also distributed under the same BSD license.
 *    Copyrights to that portion of code are retained by the U.S. government.
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
 * 3. Neither the name of the copyright holders nor the names of its contributors
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

#include "../common/common-includes.hpp"

#include "ucl/atomicutil.hpp"

#ifndef PROBLEM_SIZE
#define PROBLEM_SIZE (40)
#endif /* PROBLEM_SIZE */

// CUTOFF 0              -> no cutoff
//        1              -> elides last task-level
//        ...
//        PROBLEM_SIZE-1 -> two tasks
//        PROBLEM_SIZE   -> single task
#define CUTOFF (0)

#define PRINT_STATS 0


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


#if OMP_VERSION || SEQ_VERSION

//~ histogram<long double> hist(0.0, 1.0, 100);

fib_type partialresult;
#pragma omp threadprivate(partialresult)

template <class I>
void fib_task(I task)
{
  if (task <= 1) { partialresult += task; return; }

  while (--task > 1)
  {
    if ((CUTOFF == 0) || (task >= CUTOFF))
    {
      #pragma omp task
      fib_task<I>(task-1);
    }
    else
    {
      fib_task<I>(task-1);
    }
  }

  partialresult += 1;
}

template <class I>
auto fib_task(size_t numthreads, I num) -> std::pair<I, size_t>
{
  I                                                    res  = 0;
  size_t                                               ctr  = 1;

  #pragma omp parallel num_threads(numthreads) firstprivate(num) shared(res,ctr)
  {
    partialresult = I();

    #pragma omp single
    #pragma omp taskgroup
    {
      fib_task(num);
    }

    #pragma omp atomic
    res += partialresult;
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

template <class G, class I>
auto compute_fib(G& taskgroup, I task, ucl::simple_reducer<I>& reducer) -> void
{
  if (task <= 1) { reducer += task; return; }

  while (--task > 1)
  {
    if ((CUTOFF == 0) || (task >= CUTOFF))
      taskgroup.run( [&taskgroup, task, &reducer]()->void
                     { compute_fib(taskgroup, task-1, reducer);
                     }
                   );
    else
      compute_fib(taskgroup, task-1, reducer);
  }

  reducer += 1;
}

template <class I>
auto fib_task(size_t numthreads, I num) -> std::pair<I,size_t>
{
  tbb::task_scheduler_init init(numthreads);
  tbb::task_group          g;
  ucl::simple_reducer<I>   reducer(0);

  compute_fib(g, num, reducer);

  g.wait();
  return std::make_pair(reducer.get_value(), size_t());
}

#endif /* TBB_VERSION */


#if UCL_VERSION

struct compute_fib
{
  template <class Pool, class I>
  I operator()(Pool& pool, I task)
  {
    if (task <= 1) return task;

    I res = I(1);

    while (--task > 1)
    {
      if ((CUTOFF == 0) || (task >= CUTOFF))
        pool.enq(task-1);
      else
        res += (*this) (pool, task-1);
    }

    return res;
  }
};

template <class I>
std::pair<I, size_t> fib_task(size_t numthreads, I num)
{
  I res = ucl::execute_tasks_x(numthreads, compute_fib(), I{num});

  return std::make_pair(res, 0);
}

#endif /* UCL_VERSION */

#if CILK_VERSION

template <class I>
I compute_fib_nofork(I task)
{
  I res = I();

  while (task > 1)
  {
    res += compute_fib_nofork(task-2);

    task = task-1;
  }

  return res + task.num;
}


template <class I>
void compute_fib(I task, cilk::reducer_opadd<I>& sum)
{
  if (task <= 1) { sum += task; return; }

  while (--task > 1)
  {
    if ((CUTOFF == 0) || (task >= CUTOFF))
      cilk_spawn compute_fib(task-1, sum);
    else
    {
      compute_fib(task-1, sum);
      // faster: sum += compute_fib_nofork(fibonacci_task{task.num-2});
    }
  }

  sum += task;
}

template <class I>
std::pair<I, size_t> fib_task(size_t numthreads, I num)
{
  set_cilk_workers(numthreads);

  cilk::reducer_opadd<I> sum;

  compute_fib<I>(num, sum);
  return std::make_pair(sum.get_value(), 0);
}

#endif /* CILK_VERSION */


#if QTHREADS_VERSION

// QThreads, nowait version

template <class I>
struct qtask
{
  I          num;
  qt_sinc_t* sinc;
};

template <class I>
aligned_t compute_fib_nofork(qtask<I>& task)
{
  I         res = I();

  while (task.num > 1)
  {
    qtask<I> tsk2 = {task.num-2, task.sinc};

    res += compute_fib_nofork<I>(tsk2);
    task.num = task.num-1;
  }

  return res+task.num;
}

template <class I>
aligned_t compute_fib(void* qtsk)
{
  qtask<I>& task = *reinterpret_cast<qtask<I>*>(qtsk);

  if (task.num <= 1) { qt_sinc_submit(task.sinc, &task.num); return aligned_t(); }

  I         res = I(1);

  while (--task.num > 1)
  {
    qtask<I> tsk2 = {task.num-1, task.sinc};

    if ((CUTOFF == 0) || (task.num >= CUTOFF))
    {
      qt_sinc_expect(task.sinc, 1);
      qthread_fork_copyargs( compute_fib<I>, &tsk2, sizeof(qtask<I>), nullptr );
    }
    else
    {
      res += compute_fib_nofork<I>(tsk2);
    }
  }

  qt_sinc_submit(task.sinc, &res);
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
std::pair<I, size_t> fib_task(size_t /*numthreads*/, I num)
{
  I          result;
  qt_sinc_t* sinc   = qt_sinc_create(sizeof(I), &result, reduce<I>, 1);
  qtask<I>   task   = { num, sinc };

  qthread_fork_copyargs(compute_fib<I>, &task, sizeof(qtask<I>), nullptr);
  qt_sinc_wait(sinc, &result);
  return std::make_pair(result, 0);
}

#endif /* QTHREADS_VERSION */

#if QTHREADS_VERSION_ORIGINAL

/***
* THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION ``AS IS'' AND ANY EXPRESS OR
* IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
* MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
* EVENT SHALL SANDIA CORPORATION BE LIABLE FOR ANY DIRECT, INDIRECT,
* INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
* LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
* OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
* LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
* NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
* EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
***/

template <class I>
aligned_t fib(void *arg_)
{
    I n = *(I*)arg_;

    if (n < 2) return n;

    aligned_t ret1 = 0;
    aligned_t ret2 = 0;
    unsigned int n1 = n - 1;
    unsigned int n2 = n - 2;

    qthread_fork(fib<I>, &n1, &ret1);
    qthread_fork(fib<I>, &n2, &ret2);

    qthread_readFF(NULL, &ret1);
    qthread_readFF(NULL, &ret2);

    return ret1 + ret2;
}

/*** end SANDIA */

template <class I>
std::pair<I, size_t> fib_task(size_t /*numthreads*/, I num)
{
  aligned_t result = fib<I>(&num);

  return std::make_pair(result, 0);
}

#endif /* QTHREADS_VERSION_ORIGINAL */

int main(int argc, char** args)
{
  typedef std::chrono::time_point<std::chrono::system_clock> time_point;

  size_t   num_threads = NUMTHREADS;
  fib_type num         = PROBLEM_SIZE;

  if (argc > 1) num_threads = aux::as<size_t>(*(args+1));
  if (argc > 2) num         = aux::as<fib_type>(*(args+2));

#if QTHREADS_VERSION
  init_qthreads(num_threads);
#endif

  time_point     starttime = std::chrono::system_clock::now();

  // executes loop in parallel
  //   and uses a reduction algorithm to combine all pi values that were
  //   computed across threads.
  std::pair<fib_type, size_t> res  = fib_task(num_threads, num);

  time_point     endtime = std::chrono::system_clock::now();
  int            elapsedtime = std::chrono::duration_cast<std::chrono::milliseconds>(endtime-starttime).count();

  std::cout << "num    = " << num          << std::endl;
  std::cout << "fib    = " << res.first    << std::endl;
  std::cout << "leaves = " << res.second   << std::endl;
  std::cout << "time   = " << elapsedtime  << "ms" << std::endl;
  std::cerr << elapsedtime << std::endl;
  return 0;
}
