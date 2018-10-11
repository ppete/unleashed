/**
 * A simple implementation to solve the NQueens Problem
 *
 * (c) Peter Pirkelbauer (UAB) - 2018
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

#include "../common/common-includes.hpp"

#include "atomicutil.hpp"

#ifndef NUMTHREADS
#define NUMTHREADS (20)
#endif /* NUMTHREADS */

#ifndef PROBLEM_SIZE
#define PROBLEM_SIZE (13)
#endif /* PROBLEM_SIZE */

#define CUTOFF (PROBLEM_SIZE)
#define PRINT_STATS 0

template <size_t N>
struct board
{
  char   queens[N];
  size_t rows;

  board()
  : rows(0)
  {}
};

template <size_t N>
static inline
void print_board(board<N> task)
{
  std::cout << "N = " << N         << std::endl;
  std::cout << "r = " << task.rows << std::endl;

  for (size_t i = 0; i < task.rows; ++i)
    std::cout << size_t(task.queens[i]) << " ";

  std::cout << std::endl << std::endl;
}


template <size_t N>
static inline
bool isComplete(const board<N>& brd)
{
  return brd.rows == N;
}

template <size_t N>
static inline
bool isValid(const board<N>& board)
{
  size_t rr  = board.rows;

  if (rr < 2) return true;

  int    col = board.queens[rr-1];
  int    lhs = col;
  int    rhs = col;

  do
  {
    --rr;
    --lhs;
    ++rhs;

    int thisrow = board.queens[rr-1];

    if ((thisrow == col) || (thisrow == lhs) || (thisrow == rhs))
    {
      return false;
    }
  } while (rr > 1);

  return true;
}


#if OMP_VERSION

//~ histogram<long double> hist(0.0, 1.0, 100);

static inline
size_t thread_num()
{
  return omp_get_thread_num();
}

#if PRINT_STATS

typedef std::pair<size_t, size_t> result_t;

static inline
size_t& result(uab::aligned_type<result_t, CACHELINESZ>& v)
{
  return v.val.first;
}

static inline
size_t& count(uab::aligned_type<result_t, CACHELINESZ>& v)
{
  return v.val.second;
}

#else

typedef size_t result_t;

static inline
size_t& result(uab::aligned_type<result_t, CACHELINESZ>& v)
{
  return v.val;
}

#endif /* PRINT_STATS */

template <size_t N>
void compute_nqueens(board<N> task, uab::aligned_type<result_t, CACHELINESZ>* s)
{
  const size_t thrnum = thread_num();

  for (;;)
  {

#if PRINT_STATS
    ++(count(s[thrnum]));
#endif

    // check if last added queen is valid
    if (!isValid(task)) return;

    // if board is full
    if (isComplete(task))
    {
      ++(result(s[thrnum]));
      return;
    }

    for (size_t i = N-1; i > 0; --i)
    {
      board<N> newtask = task;

      newtask.queens[task.rows] = i;
      ++newtask.rows;

      #pragma omp task if ((CUTOFF == PROBLEM_SIZE) || (newtask.rows < CUTOFF))
      compute_nqueens(newtask, s);
    }

    task.queens[task.rows] = 0;
    ++task.rows;
  }
}

template <size_t N>
size_t nqueens_task()
{
  uab::aligned_type<result_t, CACHELINESZ> s[NUMTHREADS];
  size_t                                   res  = 0;
  board<N>                                 brd;

  omp_set_num_threads(NUMTHREADS);

  #pragma omp parallel firstprivate(brd) shared(res,s)
  {
    size_t thrnum = thread_num();

    assert(thrnum < NUMTHREADS);
    s[thrnum] = result_t();

    #pragma omp single
    #pragma omp taskgroup
    {
      compute_nqueens(brd, s);
    }

    #pragma omp atomic
    res += result(s[thrnum]);
  }

#if PRINT_STATS
  for (size_t i = 0; i < NUMTHREADS; ++i)
    std::cout << i << ": " << count(s[i]) << std::endl;
#endif

  return res;
}
#endif /* OMP_VERSION */

#if TBB_VERSION

template <class G, size_t N>
void
compute_nqueens(G& taskgroup, board<N> task, uab::simple_reducer<size_t>& reducer)
{
  for (;;)
  {
    // check if last added queen is valid
    if (!isValid(task)) return;

    // if board is full
    if (isComplete(task)) { reducer += 1; return; }

    for (size_t i = N-1; i > 0; --i)
    {
      board<N> newtask = task;

      newtask.queens[task.rows] = i;
      ++newtask.rows;

      if ((CUTOFF == PROBLEM_SIZE) || (newtask.rows < CUTOFF))
        taskgroup.run( [&taskgroup, newtask, &reducer]()->void
                     { compute_nqueens(taskgroup, newtask, reducer);
                   }
                 );
      else
        compute_nqueens(taskgroup, newtask, reducer);
    }

    task.queens[task.rows] = 0;
    ++task.rows;
  }
}

template <size_t N>
size_t nqueens_task()
{
  tbb::task_scheduler_init init(NUMTHREADS);
  tbb::task_group          g;
  uab::simple_reducer<size_t> reducer;
  board<N>                 brd;

  compute_nqueens(g, brd, reducer);

  g.wait();
  return reducer.get_value();
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


template <size_t N>
void compute_nqueens(cilk::reducer_opadd<size_t>& sum, board<N> task)
{
  for (;;)
  {
    // check if last added queen is valid
    if (!isValid(task)) return;

    // if board is full
    if (isComplete(task)) { sum += 1; return; }

    for (size_t i = N-1; i > 0; --i)
    {
      board<N> newtask = task;

      newtask.queens[task.rows] = i;
      ++newtask.rows;

      if ((CUTOFF == PROBLEM_SIZE) || (newtask.rows < CUTOFF))
      cilk_spawn compute_nqueens(sum, newtask);
      else
        compute_nqueens(sum, newtask);
    }

    task.queens[task.rows] = 0;
    ++task.rows;
  }
}

template <size_t N>
size_t nqueens_task()
{
  set_cilk_workers(NUMTHREADS);

  cilk::reducer_opadd<size_t> sum;
  board<N>                    brd;

  compute_nqueens(sum, brd);
  return sum.get_value();
}

#endif /* CILK_VERSION */


#if QTHREADS_VERSION

template <size_t N>
struct qtask
{
  board<N>    brd;
  qt_sinc_t*  sinc;
};


template <size_t N>
aligned_t compute_nqueens(void* qtsk)
{
  qtask<N>& task = *reinterpret_cast<qtask<N>*>(qtsk);

  for (;;)
  {
    // check if last added queen is valid
    if (!isValid(task.brd))
    {
      size_t res = 0;

      qt_sinc_submit(task.sinc, &res);
      delete &task;
      return aligned_t();
    }

    // if board is full
    if (isComplete(task.brd))
    {
      size_t res = 1;

      qt_sinc_submit(task.sinc, &res);
      delete &task;
      return aligned_t();
    }

    // summative increase for tasks created in loop
    qt_sinc_expect(task.sinc, N-1);

    for (size_t i = N-1; i > 0; --i)
    {
      board<N> newtask = task.brd;

      newtask.queens[task.brd.rows] = i;
      ++newtask.rows;

      if ((CUTOFF == PROBLEM_SIZE) || (newtask.rows < CUTOFF))
      qthread_fork( compute_nqueens<N>,
                    new qtask<N>{ newtask, task.sinc},
                    nullptr
                  );
      else
      {
        compute_nqueens<N>(new qtask<N>{ newtask, task.sinc});
      }
    }

    task.brd.queens[task.brd.rows] = 0;
    ++task.brd.rows;
  }
}

template <class D>
void reduce(void* target, const void* source)
{
  D*       tgt = reinterpret_cast<D*>(target);
  const D* src = reinterpret_cast<const D*>(source);

  *tgt += *src;
}

template <size_t N>
size_t nqueens_task()
{
  size_t     result = 0;
  qt_sinc_t* sinc   = qt_sinc_create(sizeof(size_t), &result, reduce<size_t>, 0);
  board<N>   brd;

  compute_nqueens<N>(new qtask<N>{brd, sinc});

  qt_sinc_wait(sinc, &result);
  return result;
}

#endif /* QTHREADS_VERSION */




#if BLAZE_VERSION

struct compute_nqueens
{
  template <class P, size_t N>
  size_t operator()(P& pool, board<N> task)
  {
    size_t res = 0;

    for (;;)
    {
      // check if last added queen is valid
      if (!isValid(task)) return res;

      // if board is full
      if (isComplete(task)) return res+1;

      for (size_t i = N-1; i > 0; --i)
      {
        board<N> newtask = task;

        newtask.queens[task.rows] = i;
        ++newtask.rows;

        if ((CUTOFF == PROBLEM_SIZE) || (newtask.rows < CUTOFF))
          pool.enq(newtask);
        else
          res += (*this)(pool, newtask);
      }

      task.queens[task.rows] = 0;
      ++task.rows;
    }
  }
};

template <size_t N>
size_t nqueens_task()
{
  board<N> brd;

  return uab::execute_tasks(NUMTHREADS, compute_nqueens(), brd);
}

#endif /* BLAZE_VERSION */

#if SEQUENTIAL_VERSION
template <size_t N>
size_t nqueens_seq(board<N> task)
{
  size_t res = 0;

  for (;;)
  {
    // check if last added queen is valid
    if (!isValid(task))
    {
      return res;
    }

    // if board is full
    if (isComplete(task)) return res+1;

    for (size_t i = N-1; i > 0; --i)
    {
      board<N> newtask = task;

      newtask.queens[task.rows] = i;
      ++newtask.rows;

      res += nqueens_seq(newtask);
    }

    task.queens[task.rows] = 0;
    ++task.rows;
  }
}

template <size_t N>
size_t nqueens_seq()
{
  return nqueens_seq(board<N>());
}
#endif /* SEQUENTIAL_VERSION */


int main()
{
#if QTHREADS_VERSION
  init_qthreads(NUMTHREADS);
#endif /* QTHREADS_VERSION */

  typedef std::chrono::time_point<std::chrono::system_clock> time_point;

  time_point     starttime = std::chrono::system_clock::now();

  // executes loop in parallel
  //   and uses a reduction algorithm to combine all pi values that were
  //   computed across threads.
  size_t num  = nqueens_task<PROBLEM_SIZE>();

  time_point     endtime = std::chrono::system_clock::now();
  int            elapsedtime = std::chrono::duration_cast<std::chrono::milliseconds>(endtime-starttime).count();

  std::cout << "num       = " << PROBLEM_SIZE << std::endl;
  std::cout << "solutions = " << num          << std::endl;
  std::cout << "time      = " << elapsedtime  << "ms" << std::endl;
  std::cerr << elapsedtime << std::endl;
  return 0;
}
