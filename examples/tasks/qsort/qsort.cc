/**
 * Task Parallel Quick Sort
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
#include <chrono>
#include <cassert>
#include <random>

#include "../common/common-includes.hpp"

#include "ucl/atomicutil.hpp"

#ifndef PROBLEM_SIZE
#define PROBLEM_SIZE 1000000
#endif /* PROBLEM_SIZE */


typedef double elem_type;

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


#if OMP_VERSION || SEQ_VERSION

template <class RandomAccessIterator, class Comparator>
void
qksort(RandomAccessIterator aa, RandomAccessIterator zz, Comparator comp)
{
  using std::swap;

  const size_t sz    = std::distance(aa, zz);

  if (sz <= 1) return;

  const size_t r = rand()%sz;

  RandomAccessIterator pivot = zz - 1;
  swap(*(aa+r), *pivot);

  RandomAccessIterator pos = aa;
  for (RandomAccessIterator it = pos; it < pivot; ++it)
  {
    if (comp(*it, *pivot))
    {
      swap(*it, *pos);
      ++pos;
    }
  }

  swap(*pos, *pivot);

  #pragma omp task firstprivate(aa, pos, comp)
  qksort(aa, pos, comp);

  qksort(pos+1, zz, comp);
}


template <class RandomAccessIterator, class Comparator>
void
quicksort(RandomAccessIterator aa, RandomAccessIterator zz, Comparator comp, size_t numthreads)
{
  #pragma omp parallel num_threads(numthreads) firstprivate(aa, zz, comp)
  #pragma omp single
  #pragma omp taskgroup
  {
    qksort(aa, zz, comp);
  }
}
#endif /* OMP_VERSION */

#if TBB_VERSION

template <class G, class RandomAccessIterator, class Comparator>
void
qksort(G& taskgroup, RandomAccessIterator aa, RandomAccessIterator zz, Comparator comp)
{
  using std::swap;

  const size_t sz    = std::distance(aa, zz);

  if (sz <= 1) return;

  const size_t r = rand()%sz;

  RandomAccessIterator pivot = zz - 1;
  swap(*(aa+r), *pivot);

  RandomAccessIterator pos = aa;
  for (RandomAccessIterator it = pos; it < pivot; ++it)
  {
    if (comp(*it, *pivot))
    {
      swap(*it, *pos);
      ++pos;
    }
  }

  swap(*pos, *pivot);

  taskgroup.run( [&taskgroup,aa,pos,comp]()->void
                 {
                   qksort(taskgroup, aa, pos, comp);
                 }
               );

  qksort(taskgroup, pos+1, zz, comp);
}


template <class RandomAccessIterator, class Comparator>
void
quicksort(RandomAccessIterator aa, RandomAccessIterator zz, Comparator comp, size_t numthreads)
{
  TBB_INIT(numthreads);
  tbb::task_group          g;

  qksort(g, aa, zz, comp);

  g.wait();
}

#endif /* TBB_VERSION */

#if CILK_VERSION

template <class RandomAccessIterator, class Comparator>
void
qksort(RandomAccessIterator aa, RandomAccessIterator zz, Comparator comp)
{
  using std::swap;

  const size_t sz    = std::distance(aa, zz);

  if (sz <= 1) return;

  const size_t r = rand()%sz;

  RandomAccessIterator pivot = zz - 1;
  swap(*(aa+r), *pivot);

  RandomAccessIterator pos = aa;
  for (RandomAccessIterator it = pos; it < pivot; ++it)
  {
    if (comp(*it, *pivot))
    {
      swap(*it, *pos);
      ++pos;
    }
  }

  swap(*pos, *pivot);

  cilk_spawn qksort(aa, pos, comp);

  qksort(pos+1, zz, comp);
}


template <class RandomAccessIterator, class Comparator>
void
quicksort(RandomAccessIterator aa, RandomAccessIterator zz, Comparator comp, size_t /*numthreads*/)
{
  //~ set_cilk_workers(numthreads);

  qksort(aa, zz, comp);
}

#endif /* CILK_VERSION */

#if QTHREADS_VERSION

#error "N/A"

#endif /* QTHREADS_VERSION */



#if UCL_VERSION

template <class RandomAccessIterator>
struct qs_task
{
  RandomAccessIterator aa;
  RandomAccessIterator zz;
};

template <class Comparator>
struct qksort
{
  Comparator comp;

  explicit
  qksort(Comparator c)
  : comp(c)
  {}

  template <class pool_t, class RandomAccessIterator>
  ucl::Void
  operator()(pool_t& pool, qs_task<RandomAccessIterator> task)
  {
    typedef qs_task<RandomAccessIterator> qs_task;

    using std::swap;

    const size_t sz    = std::distance(task.aa, task.zz);

    if (sz <= 1) return ucl::Void();

    const size_t r = rand()%sz;

    RandomAccessIterator pivot = task.zz - 1;
    swap(*(task.aa+r), *pivot);

    RandomAccessIterator pos = task.aa;
    for (RandomAccessIterator it = pos; it < pivot; ++it)
    {
      if (comp(*it, *pivot))
      {
        swap(*it, *pos);
        ++pos;
      }
    }

    swap(*pos, *pivot);

    //~ if (sz > 10000)
      pool.enq(qs_task { task.aa,    pos } );
    //~ else
      //~ (*this)( pool, qs_task { task.aa, pos } );

    return (*this)( pool, qs_task { pos+1, task.zz  } );
  }
};

template <class RandomAccessIterator, class Comparator>
void quicksort(RandomAccessIterator aa, RandomAccessIterator zz, Comparator comp, size_t numthreads)
{
  typedef qs_task<RandomAccessIterator> qs_task;

  qksort<Comparator> qs(comp);

  ucl::execute_tasks_x(numthreads, qs, qs_task{ aa, zz });
}

#endif /* UCL_VERSION */

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
  quicksort(problem, problem+problem_sz, std::less<elem_type>(), num_threads);

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
