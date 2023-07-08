/**
 * A simple implementation to solve the NQueens Problem
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
#include <thread>
#include <list>
#include <limits>

#include "../common/common-includes.hpp"

  //~ #include "ucl/atomicutil.hpp"

#ifndef NUMTHREADS
#define NUMTHREADS (20)
#endif /* NUMTHREADS */

#ifndef PROBLEM_SIZE
#define PROBLEM_SIZE (13)
#endif /* PROBLEM_SIZE */

std::size_t global_cut_off = 0;


// we use a dynamic_char_seq instead of a std::string to avoid
//   uncertainties with implementation details.
//   e.g., a copy does not need to have the same capacity as the original string,
//         small string optimizations, ...
struct dynamic_char_seq
{
    explicit
    dynamic_char_seq(std::size_t max_capacity)
    : capa(max_capacity), count(0), data(new char[capa])
    {}

    dynamic_char_seq(const dynamic_char_seq& orig)
    : capa(orig.capacity()), count(orig.size()), data(new char[capa])
    {
      std::copy(&orig.data[0], &orig.data[orig.size()], &data[0]);
    }

    dynamic_char_seq& operator=(const dynamic_char_seq& orig)
    {
      dynamic_char_seq cpy(orig);

      *this = std::move(cpy);
      return *this;
    }

    dynamic_char_seq(dynamic_char_seq&& orig)            = default;
    dynamic_char_seq& operator=(dynamic_char_seq&& orig) = default;

    std::size_t size() const { return count; }
    std::size_t capacity() const { return capa; }

    std::size_t at(std::size_t pos) const
    {
      assert(pos < count);

      return data[pos];
    }

    void push_back(char ch)
    {
      assert(count < capa);
      data[count] = ch;
      ++count;
    }

  private:
    std::size_t             capa;
    std::size_t             count;
    std::unique_ptr<char[]> data;
};

struct board
{
    explicit
    board(size_t sz)
    : queens(sz)
    {}

    board(const board&)            = default;
    board(board&&)                 = default;

    board& operator=(const board&) = default;
    board& operator=(board&&)      = default;

    size_t size()   const            { return queens.capacity(); }
    size_t rows()   const            { return queens.size(); }
    bool   complete() const          { return queens.capacity() == queens.size(); }
    size_t at(std::size_t row) const { return queens.at(row); }

    void   append(char pos)          { queens.push_back(pos); }

    bool   valid()    const;

  private:
    dynamic_char_seq queens;

    board() = delete;
};


void print_board(const board& task)
{
  std::cout << "N = " << task.size() << std::endl;
  std::cout << "r = " << task.rows() << std::endl;

  for (size_t i = 0; i < task.rows(); ++i)
    std::cout << task.at(i) << " ";

  std::cout << std::endl << std::endl;
}


bool board::valid() const
{
  size_t rr = rows();

  if (rr < 2) return true;

  int    col = queens.at(rr-1);
  int    lhs = col;
  int    rhs = col;

  do
  {
    --rr;
    --lhs;
    ++rhs;

    int thisrow = queens.at(rr-1);

    if ((thisrow == col) || (thisrow == lhs) || (thisrow == rhs))
    {
      return false;
    }
  } while (rr > 1);

  return true;
}


#if OMP_VERSION || SEQ_VERSION

size_t partialresult;
#pragma omp threadprivate(partialresult)

void compute_nqueens(board task)
{
  for (;;)
  {
    // check if last added queen is valid
    if (!task.valid()) return;

    // if board is full
    if (task.complete()) { ++partialresult; return; }

    for (size_t i = task.size()-1; i > 0; --i)
    {
      board newtask = task;

      newtask.append(i);

      #pragma omp task if (newtask.rows() <= global_cut_off)
      compute_nqueens(newtask);
    }

    task.append(0);
  }
}


size_t nqueens_task(size_t numthreads, size_t problem_size)
{
  // avoid unused parameter warning for sequential version
  size_t res = (numthreads-numthreads);

  #pragma omp parallel num_threads(numthreads) shared(res)
  {
    #pragma omp single
    #pragma omp taskgroup
    {
      compute_nqueens(board{problem_size});
    }

    #pragma omp atomic
    res += partialresult;
  }

  return res;
}
#endif /* OMP_VERSION */

#if TBB_VERSION

template <class G>
void
compute_nqueens(G& taskgroup, const board& task, ucl::simple_reducer<size_t>& reducer)
{
  // check if last added queen is valid
  if (!task.valid()) return;

  // if board is full
  if (task.complete()) { reducer += 1; return; }

  for (int i = task.size(); i >= 0; --i)
  {
    board newtask = task;

    newtask.append(i);

    if (newtask.rows() <= global_cut_off)
      taskgroup.run( [&taskgroup, t = std::move(newtask), &reducer]()->void
                     {
										   compute_nqueens(taskgroup, t, reducer);
                     }
                   );
    else
      compute_nqueens(taskgroup, newtask, reducer);
  }
}

size_t nqueens_task(size_t numthreads, size_t problem_size)
{
  TBB_INIT(numthreads);
  TBB::task_group             g;
  ucl::simple_reducer<size_t> reducer;

  compute_nqueens(g, board{problem_size}, reducer);

  g.wait();
  return reducer.get_value();
}

#endif /* TBB_VERSION */

#if CILK_VERSION

void count_init(void* sum) { *static_cast<std::size_t*>(sum) = 0.0; }
void count_plus(void* lhs, void* rhs) { *static_cast<std::size_t*>(lhs) += *static_cast<std::size_t*>(rhs); }

std::size_t cilk_reducer(count_init, count_plus) count(0);

void compute_nqueens(board task)
{
  for (;;)
  {
    // check if last added queen is valid
    if (!task.valid()) return;

    // if board is full
    if (task.complete()) { count += 1; return; }

    for (size_t i = task.size()-1; i > 0; --i)
    {
      board newtask = task;

      newtask.append(i);

      if (newtask.rows() <= global_cut_off)
        cilk_spawn compute_nqueens(newtask);
      else
        compute_nqueens(newtask);
    }

    task.append(0);
  }
}

size_t nqueens_task(size_t /*numthreads*/, size_t problem_size)
{
  //~ set_cilk_workers(numthreads);

  //~ cilk::reducer_opadd<size_t> sum;

  compute_nqueens(board{problem_size});
  return count;
}

#endif /* CILK_VERSION */



#if UCL_VERSION

struct compute_nqueens
{
  template <class P>
  size_t operator()(P& pool, board task)
  {
    size_t res = 0;

    for (;;)
    {
      // check if last added queen is valid
      if (!task.valid()) return res;

      // if board is full
      if (task.complete()) return res+1;

      for (size_t i = task.size()-1; i > 0; --i)
      {
        board newtask = task;

        newtask.append(i);

        if (newtask.rows() <= global_cut_off)
          pool.enq(newtask);
        else
          res += (*this)(pool, newtask);
      }

      // local task
      task.append(0);
    }
  }
};

size_t nqueens_task(size_t numthreads, size_t problem_size)
{
  return ucl::execute_tasks_x(numthreads, compute_nqueens{}, board{problem_size});
}

#endif /* UCL_VERSION */


int main(int argc, char** args)
{
  typedef std::chrono::time_point<std::chrono::system_clock> time_point;

  size_t   num_threads  = NUMTHREADS;
  size_t   problem_size = PROBLEM_SIZE;

  if (argc > 1) num_threads  = aux::as<size_t>(*(args+1));
  if (argc > 2) problem_size = aux::as<size_t>(*(args+2));

  global_cut_off = (argc > 3) ? aux::as<size_t>(*(args+3)) : problem_size;

  if (problem_size > std::numeric_limits<char>::max())
  {
    std::cerr << "problem size must be <= " << std::numeric_limits<char>::max()
              << std::endl;
    exit(0);
  }

#if QTHREADS_VERSION
  init_qthreads(num_threads);
#endif

  time_point     starttime = std::chrono::system_clock::now();

  const size_t   num  = nqueens_task(num_threads, problem_size);

  time_point     endtime = std::chrono::system_clock::now();
  int            elapsedtime = std::chrono::duration_cast<std::chrono::milliseconds>(endtime-starttime).count();

  std::cout << "threads      = " << num_threads << std::endl;
  std::cout << "problem size = " << problem_size << std::endl;
  std::cout << "cutoff       = " << global_cut_off << (problem_size <= global_cut_off ? " (no cutoff)" : "") << std::endl;
  std::cout << "solutions    = " << num          << std::endl;
  std::cout << "time         = " << elapsedtime  << "ms" << std::endl;
  std::cerr << elapsedtime << std::endl;
  return 0;
}
