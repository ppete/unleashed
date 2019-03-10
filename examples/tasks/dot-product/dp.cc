/// Task-based approach for computing A . B / |A| * |B|
///
/// Implementer: Peter Pirkelbauer

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
#include <cassert>
#include <random>

#include "../common/common-includes.hpp"

#include "ucl/atomicutil.hpp"

// Problem size
#ifndef PROBLEM_SIZE
#define PROBLEM_SIZE (100000000)
#endif /* PROBLEM_SIZE */

typedef long double product_type;

template <class D>
struct dp_result
{
  dp_result()
  : product(0), lhs_len_sq(0), rhs_len_sq(0)
  {}

  dp_result(D prod, D llensq, D rlensq)
  : product(prod), lhs_len_sq(llensq), rhs_len_sq(rlensq)
  {}

  dp_result& operator+=(const dp_result& rhs)
  {
    this->product    += rhs.product;
    this->lhs_len_sq += rhs.lhs_len_sq;
    this->rhs_len_sq += rhs.rhs_len_sq;

    return *this;
  }

  D   product;
  D   lhs_len_sq;
  D   rhs_len_sq;
};


#if TBB_VERSION

template <class G, class D>
void dp_compute( G& taskgroup,
                 D* lhs,
                 D* rhs,
                 size_t lo,
                 size_t hi,
                 ucl::simple_reducer<dp_result<D> >& rdcr
               )
{
  while (lo + 1 < hi)
  {
    size_t mid = (hi+lo) / 2;

    if (mid < hi)
    {
      taskgroup.run( [&taskgroup, lhs, rhs, mid, hi, &rdcr]()->void
                     {
                       dp_compute(taskgroup, lhs, rhs, mid, hi, rdcr);
                     }
                   );
    }

    hi = mid;
  }

  rdcr += dp_result<D>{ lhs[lo] * rhs[lo], lhs[lo] * lhs[lo], rhs[lo] * rhs[lo] };
}


template <class D>
D dp_calc(D* lhs, D* rhs, size_t numthreads, size_t len)
{
  tbb::task_scheduler_init           init(numthreads);
  tbb::task_group                    g;
  ucl::simple_reducer<dp_result<D> > reducer;

  dp_compute(g, lhs, rhs, 0, len, reducer);
  g.wait();

  dp_result<D>                       res = reducer.get_value();

  return res.product / (std::sqrt(res.lhs_len_sq) * std::sqrt(res.rhs_len_sq));
}

#endif /* TBB_VERSION */

#if UCL_VERSION

template <class D>
struct dp_task
{
  size_t lo;
  size_t hi;
};

template <class D>
struct dp_operator
{
  D*     lhs;
  D*     rhs;

  template <class P>
  auto operator()(P& pool, dp_task<D> task) -> dp_result<D>
  {
    while (task.lo + 1 < task.hi)
    {
      size_t mid = (task.hi+task.lo) / 2;

      if (mid < task.hi)
        pool.enq(dp_task<D>{ mid, task.hi } );

      task.hi = mid;
    }

    return dp_result<D> { lhs[task.lo] * rhs[task.lo],
                          lhs[task.lo] * lhs[task.lo],
                          rhs[task.lo] * rhs[task.lo]
                        };
  }
};

template <class D>
auto dp_calc(D* lhs, D* rhs, size_t numthreads, size_t len) -> D
{
  dp_result<D> res = ucl::execute_tasks_x(numthreads, dp_operator<D>{lhs,rhs}, dp_task<D>{ 0, len });

  return res.product / (std::sqrt(res.lhs_len_sq) * std::sqrt(res.rhs_len_sq));
}
#endif /* UCL_VERSION */

#if OMP_VERSION

typedef dp_result<product_type> reduction_type;

static product_type red_product;
static product_type red_lhs_sum_sq;
static product_type red_rhs_sum_sq;

#pragma omp threadprivate(red_product, red_lhs_sum_sq, red_rhs_sum_sq)

template <class D>
void dp_compute(D* lhs, D* rhs, size_t lo, size_t hi)
{
  while (lo + 1 < hi)
  {
    size_t mid = (hi+lo) / 2;

    if (mid < hi)
    {
      #pragma omp task
      dp_compute(lhs, rhs, mid, hi);
    }

    hi = mid;
  }

  red_product    += lhs[lo] * rhs[lo];
  red_lhs_sum_sq += lhs[lo] * lhs[lo];
  red_rhs_sum_sq += rhs[lo] * rhs[lo];
}


template <class D>
auto dp_calc(D* lhs, D* rhs, size_t numthreads, size_t len) -> D
{
  reduction_type res;

  #pragma omp parallel num_threads(numthreads) firstprivate(lhs, rhs, len)
  {
    red_product = red_lhs_sum_sq = red_rhs_sum_sq = 0;

    #pragma omp single
    #pragma omp taskgroup
    {
      dp_compute(lhs, rhs, 0, len);
    }

    #pragma omp atomic
    res.product    += red_product;

    #pragma omp atomic
    res.lhs_len_sq += red_lhs_sum_sq;

    #pragma omp atomic
    res.rhs_len_sq += red_rhs_sum_sq;
  }

  return res.product / (std::sqrt(res.lhs_len_sq) * std::sqrt(res.rhs_len_sq));
}
#endif /* OMP_VERSION */

#if CILK_VERSION

template <class D>
void dp_compute( D* lhs, D* rhs, size_t lo, size_t hi,
                 cilk::reducer_opadd<D>& product,
                 cilk::reducer_opadd<D>& lhs_len_sq,
                 cilk::reducer_opadd<D>& rhs_len_sq
               )
{
  while (lo + 1 < hi)
  {
    size_t mid = (hi+lo) / 2;

    if (mid < hi)
    {
      cilk_spawn dp_compute<D>(lhs, rhs, mid, hi, product, lhs_len_sq, rhs_len_sq);
    }

    hi = mid;
  }

  product     += lhs[lo] * rhs[lo];
  lhs_len_sq  += lhs[lo] * lhs[lo];
  rhs_len_sq  += rhs[lo] * rhs[lo];
}


template <class D>
auto dp_calc(D* lhs, D* rhs, size_t numthreads, size_t len) -> D
{
  set_cilk_workers(numthreads);

  cilk::reducer_opadd<D> product;
  cilk::reducer_opadd<D> lhs_len_sq;
  cilk::reducer_opadd<D> rhs_len_sq;

  dp_compute(lhs, rhs, 0, len, product, lhs_len_sq, rhs_len_sq);

  return product.get_value() / (std::sqrt(lhs_len_sq.get_value()) * std::sqrt(rhs_len_sq.get_value()));
}

#endif /* CILK_VERSION */


#if QTHREADS_VERSION

template <class D>
struct dp_task
{
  size_t     lo;
  size_t     hi;
  D*         lhs;
  D*         rhs;
  qt_sinc_t* sinc;
};


template <class D>
aligned_t dp_compute(void* qtsk)
{
  typedef dp_task<D> dp_task;

  dp_task& task = *reinterpret_cast<dp_task*>(qtsk);

  while (task.lo + 1 < task.hi)
  {
    size_t mid = (task.hi+task.lo) / 2;

    if (mid < task.hi)
    {
      qt_sinc_expect(task.sinc, 1);
      qthread_fork( dp_compute<D>,
                    new dp_task{ mid, task.hi, task.lhs, task.rhs, task.sinc},
                    nullptr
                  );

    }

    task.hi = mid;
  }

  dp_result<D> res = { task.lhs[task.lo] * task.rhs[task.lo],
                       task.lhs[task.lo] * task.lhs[task.lo],
                       task.rhs[task.lo] * task.rhs[task.lo]
                     };
  qt_sinc_submit(task.sinc, &res);
  delete &task;
  return 0;
}

template <class D>
void reduce(void* target, const void* source)
{
  dp_result<D>*       tgt = reinterpret_cast<dp_result<D>*>(target);
  const dp_result<D>* src = reinterpret_cast<const dp_result<D>*>(source);

  *tgt += *src;
}

template <class D>
auto dp_calc(D* lhs, D* rhs, size_t /*numthreads*/, size_t len) -> D
{
  dp_result<D> res { 0, 0, 0 };
  qt_sinc_t*   sinc = qt_sinc_create(sizeof(dp_result<D>), &res, reduce<D>, 1);

  dp_compute<D>(new dp_task<D>{ 0, len, lhs, rhs, sinc});
  qt_sinc_wait(sinc, &res);

  return res.product / (std::sqrt(res.lhs_len_sq) * std::sqrt(res.rhs_len_sq));
}

#endif /* QTHREADS_VERSION */

template <class D>
void init_vectors(D* lhs, D* rhs, size_t len)
{
  std::uniform_real_distribution<D> dist( -100, 100 );
  std::default_random_engine        rand(3);

  while (len--)
  {
    lhs[len] = dist(rand);
    rhs[len] = dist(rand);
  }
}

int main(int argc, char** args)
{
  size_t num_threads = NUMTHREADS;
  size_t num_elems   = PROBLEM_SIZE;

  if (argc > 1) num_threads = aux::as<size_t>(*(args+1));
  if (argc > 2) num_elems   = aux::as<size_t>(*(args+2));


#if QTHREADS_VERSION
  init_qthreads(NUMTHREADS);
#endif /* QTHREADS_VERSION */

  product_type* lhs = new product_type[num_elems];
  product_type* rhs = new product_type[num_elems];

  init_vectors(lhs, rhs, num_elems);

  typedef std::chrono::time_point<std::chrono::system_clock> time_point;
  time_point     starttime = std::chrono::system_clock::now();

  product_type res = dp_calc(lhs, rhs, num_threads, num_elems);

  time_point     endtime = std::chrono::system_clock::now();
  int            elapsedtime = std::chrono::duration_cast<std::chrono::milliseconds>(endtime-starttime).count();

  std::cout << "time = " << elapsedtime << "ms" << std::endl;
  std::cout << "Final result = " << res << std::endl;
  std::cerr << elapsedtime << std::endl;

  delete[] lhs;
  delete[] rhs;
  return 0;
}
