/**
 * Implements merge sort
 *
 * Implementer: Peter Pirkelbauer (LLNL)
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
#include <chrono>
#include <cassert>
#include <random>
#include <vector>

#include "../common/common-includes.hpp"

#if STL_VERSION

#include <execution>
#include <algorithm>

#else

#include "ucl/atomicutil.hpp"

#endif /* STL_VERSION */

#ifndef PROBLEM_SIZE
#define PROBLEM_SIZE 1000000
#endif /* PROBLEM_SIZE */

typedef long double elem_type;

template <class RandomAccessIterator, class Comparator>
bool validate_array(RandomAccessIterator aa, RandomAccessIterator zz, Comparator comp)
{
  if (aa == zz) return true;

  RandomAccessIterator prv = aa;

  ++aa;
  while (aa < zz && comp(*prv, *aa))
  {
    prv = aa;
    ++aa;
  }

  return aa == zz;
}


template <class RandomAccessIterator, class Comparator>
void merge0( RandomAccessIterator aa,
             RandomAccessIterator mi,
             RandomAccessIterator zz,
             Comparator comp
           )
{
  typedef typename std::iterator_traits<RandomAccessIterator>::value_type value_type;


  const size_t               num   = std::distance(aa, zz);
  const RandomAccessIterator start = aa;
  const RandomAccessIterator md    = mi;
  std::vector<value_type>    tmp;

  tmp.reserve(num);
  while (aa < md && mi < zz)
  {
    if (comp(*aa, *mi))
    {
      tmp.emplace_back(std::move(*aa));
      ++aa;
    }
    else
    {
      tmp.emplace_back(std::move(*mi));
      ++mi;
    }
  }

  while (aa < md)
  {
    tmp.emplace_back(std::move(*aa));
    ++aa;
  }

  while (mi < zz)
  {
    tmp.emplace_back(std::move(*mi));
    ++mi;
  }

  std::move(tmp.begin(), tmp.end(), start);
}


#if WOMP_VERSION || SEQ_VERSION

template <class RandomAccessIterator, class Comparator>
void mgsort(RandomAccessIterator aa, RandomAccessIterator zz, Comparator comp)
{
  const size_t         diff = std::distance(aa, zz);

  if (diff <= 1) return;

  RandomAccessIterator mid  = aa+diff/2;

  #pragma omp task firstprivate(aa, mid, comp)
  mgsort(aa, mid, comp);

  mgsort(mid, zz, comp);

  #pragma omp taskwait

  merge0(aa, mid, zz, comp);
}

template <class RandomAccessIterator, class Comparator>
void
mergesort(RandomAccessIterator aa, RandomAccessIterator zz, Comparator comp, size_t numthreads)
{
  #pragma omp parallel num_threads(numthreads) firstprivate(aa, zz, comp)
  #pragma omp single
  #pragma omp taskgroup
  {
    mgsort(aa, zz, comp);
  }
}

#endif /* WOMP_VERSION */

#if OMP_VERSION

typedef ucl::continuation task_continuation;
// typedef omp_continue task_continuation;

template <class RandomAccessIterator, class Comparator>
void mgsort(RandomAccessIterator aa, RandomAccessIterator zz, Comparator comp, task_continuation parent)
{
  const size_t         diff = std::distance(aa, zz);

  if (diff <= 1) return;

  RandomAccessIterator mid  = aa+diff/2;
  task_continuation    cont( [ aa, mid, zz, comp, parent]() -> void
                             {
                               merge0(aa, mid, zz, comp);
                             },
                             1
                           );

  #pragma omp task firstprivate(aa, mid, cont)
  mgsort(aa, mid, comp, cont);

  mgsort(mid, zz, comp, cont);
}

template <class RandomAccessIterator, class Comparator>
void
mergesort(RandomAccessIterator aa, RandomAccessIterator zz, Comparator comp, size_t numthreads)
{
  #pragma omp parallel num_threads(numthreads) firstprivate(aa, zz, comp)
  #pragma omp single
  #pragma omp taskgroup
  {
    task_continuation nil;

    mgsort(aa, zz, comp, std::move(nil));
  }
}

#endif /* OMP_VERSION */

#if TBB_VERSION

#endif /* TBB_VERSION */

#if CILK_VERSION

template <class RandomAccessIterator, class Comparator>
void mgsort(RandomAccessIterator aa, RandomAccessIterator zz, Comparator comp)
{
  const size_t         diff = std::distance(aa, zz);

  if (diff <= 1) return;

  RandomAccessIterator mid  = aa+diff/2;

  cilk_spawn mgsort(aa, mid, comp);

  mgsort(mid, zz, comp);

  cilk_sync;

  merge0(aa, mid, zz, comp);
}


template <class RandomAccessIterator, class Comparator>
void
mergesort(RandomAccessIterator aa, RandomAccessIterator zz, Comparator comp, size_t numthreads)
{
  set_cilk_workers(numthreads);

  mgsort(aa, zz, comp);
}


#endif /* CILK_VERSION */

#if QTHREADS_VERSION

#error "N/A"

#endif /* QTHREADS_VERSION */



#if UCL_VERSION

template <class RandomAccessIterator>
struct ms_task
{
  RandomAccessIterator aa;
  RandomAccessIterator zz;
  ucl::continuation    cont;

  ms_task()
  : aa(), zz(), cont()
  {}

  ms_task(RandomAccessIterator start, RandomAccessIterator limit, ucl::continuation&& c)
  : aa(start), zz(limit), cont(std::move(c))
  {}
};

template <class Comparator>
struct mgsort
{
  Comparator comp;

  explicit
  mgsort(Comparator c)
  : comp(c)
  {}

  template <class pool_t, class RandomAccessIterator>
  ucl::Void
  operator()(pool_t& pool, ms_task<RandomAccessIterator> t)
  {
    typedef ms_task<RandomAccessIterator> ms_task;

    RandomAccessIterator aa   = t.aa;
    RandomAccessIterator zz   = t.zz;
    const size_t         diff = std::distance(aa, zz);

    if (diff <= 1) return ucl::Void();

    RandomAccessIterator mid  = aa+diff/2;
    ucl::continuation    cont( [ aa, mid, zz, cmpfn = this->comp, cont = std::move(t.cont) ]() -> void
                               {
                                 merge0(aa, mid, zz, cmpfn);
                               },
                               2
                             );

    pool.enq(ms_task{aa, mid, cont.copy_external_count()});
    return (*this)(pool, ms_task{mid, zz, std::move(cont)});
  }
};

template <class RandomAccessIterator, class Comparator>
void mergesort(RandomAccessIterator aa, RandomAccessIterator zz, Comparator comp, size_t numthreads)
{
  typedef ms_task<RandomAccessIterator> ms_task;

  mgsort<Comparator> ms(comp);
  ucl::continuation  nil;
  ms_task            task(aa, zz, std::move(nil));

  ucl::execute_tasks_x(numthreads, ms, std::move(task));
}

#endif /* UCL_VERSION */

#if STL_VERSION

template <class RandomAccessIterator, class Comparator>
void mergesort(RandomAccessIterator aa, RandomAccessIterator zz, Comparator comp, size_t /*numthreads*/)
{
  std::sort(std::execution::par_unseq, aa, zz, comp);
}

#endif /* STL_VERSION */



template <class ForwardIterator>
void init_array(ForwardIterator aa, ForwardIterator zz)
{
  typedef typename std::iterator_traits<ForwardIterator>::value_type value_type;

  std::uniform_real_distribution<value_type> dist( -100, 100 );
  std::default_random_engine                 rand(3);

  while (aa != zz)
  {
    *aa = dist(rand);
    ++aa;
  }
}

int main(int argc, char** args)
{
  typedef std::chrono::time_point<std::chrono::system_clock> time_point;

  size_t num_threads = NUMTHREADS;
  size_t problem_sz  = PROBLEM_SIZE;

  if (argc > 1) num_threads = aux::as<size_t>(*(args+1));
  if (argc > 2) problem_sz  = aux::as<size_t>(*(args+2));

#if QTHREADS_VERSION
  init_qthreads(num_threads);
#endif

  elem_type* problem = new elem_type[problem_sz];

  init_array(problem, problem+problem_sz);

  time_point     starttime = std::chrono::system_clock::now();

  // executes task-based parallel qksort
  mergesort(problem, problem+problem_sz, std::less<elem_type>(), num_threads);

  time_point     endtime = std::chrono::system_clock::now();
  int            elapsedtime = std::chrono::duration_cast<std::chrono::milliseconds>(endtime-starttime).count();

  if (!validate_array(problem, problem+problem_sz, std::less<elem_type>()))
  {
    std::cout << "oops, not sorted" << std::endl;
    std::cerr << "oops, not sorted" << std::endl;
  }
  else
  {
    std::cout << "time = " << elapsedtime << "ms" << std::endl;
    std::cerr << elapsedtime << std::endl;
  }

  delete[] problem;
  return 0;
}
