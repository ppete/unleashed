/// Task-based solution for the Traveling Salesman Problem
///
/// Implementer: Christina Peterson (UCF), Peter Pirkelbauer (UAB)

/**
 * This program is part of the Blaze-Task Test Suite
 * Copyright (c) 2018, University of Central Florida
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


#include <stdio.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <memory>

#ifndef NUMTHREADS
#define NUMTHREADS (20)
#endif /* NUMTHREADS */

// Problem size
#ifndef PROBLEM_SIZE
#define PROBLEM_SIZE (11)
#endif /* PROBLEM_SIZE */

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

#define MAX_NODES 12

typedef long result_t;

std::atomic<result_t> min_path;

struct BranchSet
{
  int left_branch[PROBLEM_SIZE][PROBLEM_SIZE];
  int right_branch[PROBLEM_SIZE][PROBLEM_SIZE];
};

struct tsp_task
{
  typedef std::shared_ptr<BranchSet> branch_ptr;
  //~ typedef BranchSet* branch_ptr;

  int**      graph;
  branch_ptr branch;
  bool       left;

  tsp_task()
  : graph(nullptr), branch(nullptr), left(false)
  {}

  explicit
  tsp_task(int** g, branch_ptr b, bool l)
  : graph(g), branch(b), left(l)
  {}
};

void print_graph(int graph[PROBLEM_SIZE][PROBLEM_SIZE])
{
  for(size_t i = 0; i < PROBLEM_SIZE; i++)
  {
    for(size_t j = 0; j < PROBLEM_SIZE; j++)
    {
      printf("%d ", graph[i][j]);
    }
    printf("\n");
  }
}

enum branch_kind
{
  INCLUDE = 1,
  EXCLUDE   = 2,
  UNDECIDED = 3
};


static inline
void upd_lower_bound_if_needed(result_t val)
{
  result_t min = min_path.load(std::memory_order_relaxed);

  while (val < min)
  {
    min_path.compare_exchange_strong(min, val, std::memory_order_relaxed, std::memory_order_relaxed);
  }
}


void left_include_edges(BranchSet* branch, int* include_count, int* exclude_count, size_t i);
void left_exclude_edges(BranchSet* branch, int* include_count, int* exclude_count, size_t i);

void left_include_edges(BranchSet* branch, int* include_count, int* exclude_count, size_t i)
{
  //printf("Inside left_include_edges\n");
  for(size_t j = 0; j < PROBLEM_SIZE; ++j)
  {
    if(j != i)
    {
      if(branch->left_branch[i][j] == UNDECIDED)
      {
        branch->left_branch[i][j] = INCLUDE;
        branch->left_branch[j][i] = INCLUDE;

        include_count[i] = include_count[i] + 1;
        include_count[j] = include_count[j] + 1;

        if(include_count[j] == 2)
        {
          left_exclude_edges(branch, include_count, exclude_count, j);
        }
      }
    }
  }
}

void left_exclude_edges(BranchSet* branch, int* include_count, int* exclude_count, size_t i)
{
  //printf("Inside left_exclude_edges\n");
  for(size_t j = 0; j < PROBLEM_SIZE; ++j)
  {
    if(j != i)
    {
      if(branch->left_branch[i][j] == UNDECIDED)
      {
        branch->left_branch[i][j] = EXCLUDE;
        branch->left_branch[j][i] = EXCLUDE;

        exclude_count[i] = exclude_count[i] + 1;
        exclude_count[j] = exclude_count[j] + 1;

        if(exclude_count[j] == (PROBLEM_SIZE - 3))
        {
          left_include_edges(branch, include_count, exclude_count, j);
        }
      }
    }
  }
}

void right_include_edges(BranchSet* branch, int* include_count, int* exclude_count, size_t i);
void right_exclude_edges(BranchSet* branch, int* include_count, int* exclude_count, size_t i);

void right_include_edges(BranchSet* branch, int* include_count, int* exclude_count, size_t i)
{
  //printf("Inside right_include_edges\n");
  for (size_t j = 0; j < PROBLEM_SIZE; ++j)
  {
    if(j != i)
    {
      if(branch->right_branch[i][j] == UNDECIDED)
      {
        branch->right_branch[i][j] = INCLUDE;
        branch->right_branch[j][i] = INCLUDE;

        include_count[i] = include_count[i] + 1;
        include_count[j] = include_count[j] + 1;

        if(include_count[j] == 2)
        {
          right_exclude_edges(branch, include_count, exclude_count, j);
        }
      }
    }
  }
}

void right_exclude_edges(BranchSet* branch, int* include_count, int* exclude_count, size_t i)
{
  //printf("Inside right_exclude_edges\n");
  for (size_t j = 0; j < PROBLEM_SIZE; ++j)
  {
    if(j != i)
    {
      if(branch->right_branch[i][j] == UNDECIDED)
      {
        branch->right_branch[i][j] = EXCLUDE;
        branch->right_branch[j][i] = EXCLUDE;

        exclude_count[i] = exclude_count[i] + 1;
        exclude_count[j] = exclude_count[j] + 1;

        if(exclude_count[j] == (PROBLEM_SIZE - 3))
        {
          right_include_edges(branch, include_count, exclude_count, j);
        }
      }
    }
  }
}


BranchSet* create_branch(int explore[PROBLEM_SIZE][PROBLEM_SIZE])
{
  BranchSet* branches = new BranchSet;

  int left_inc_count[PROBLEM_SIZE];
  int left_exc_count[PROBLEM_SIZE];

  int right_inc_count[PROBLEM_SIZE];
  int right_exc_count[PROBLEM_SIZE];

  //Try to find a more efficient way of handling this later
  for (size_t i = 0; i < PROBLEM_SIZE; ++i)
  {
    left_inc_count[i] = 0;
    left_exc_count[i] = 0;
    right_inc_count[i] = 0;
    right_exc_count[i] = 0;

    for (size_t j = 0; j < PROBLEM_SIZE; ++j)
    {
      branches->left_branch[i][j] = explore[i][j];
      branches->right_branch[i][j] = explore[i][j];

      if(explore[i][j] == INCLUDE)
      {
        left_inc_count[i] = left_inc_count[i] + 1;
        right_inc_count[i] = right_inc_count[i] + 1;
      } else if (explore[i][j] == EXCLUDE) {
        left_exc_count[i] = left_exc_count[i] + 1;
        right_exc_count[i] = right_exc_count[i] + 1;
      } else {

      }
    }
  }

  bool done = false;

  for(size_t i = 0; i < PROBLEM_SIZE; ++i)
  {
    for(size_t j = 0; j < PROBLEM_SIZE; ++j)
    {
      if(i == j)
      {
        branches->left_branch[i][i] = UNDECIDED;
        branches->right_branch[i][i] = UNDECIDED;
      }
      else
      {
        // A decision has not yet been made regarding this edge
        if (explore[i][j] == UNDECIDED && done == false)
        {
          //Set left_branch
          branches->left_branch[i][j] = INCLUDE; //This edge will be included in the explored path
          branches->left_branch[j][i] = INCLUDE;

          //printf("Setting branches->left_branch[%lu][%lu] = %d\n", i, j, branches->left_branch[j][i]);
          left_inc_count[i] = left_inc_count[i] + 1;
          left_inc_count[j] = left_inc_count[j] + 1;

          if(left_inc_count[i] == 2)
          {
            //printf("Calling left_exclude_edges\n");
            left_exclude_edges(branches, left_inc_count, left_exc_count, i);
          }

          if(left_inc_count[j] == 2)
          {
            //printf("Calling left_exclude_edges\n");
            left_exclude_edges(branches, left_inc_count, left_exc_count, j);
          }

          if(left_inc_count[i] > 2 || left_inc_count[j] > 2)
          {
            //printf("Error: left_inc_count > 2\n");
            delete branches;
            return nullptr;
          }

          //Set right_branch
          branches->right_branch[i][j] = EXCLUDE; //This edge will not be included in the explored path
          branches->right_branch[j][i] = EXCLUDE;
          //printf("Setting branches->right_branch[%lu][%lu] = %d\n", i, j, branches->right_branch[j][i]);
          right_exc_count[i] = right_exc_count[i] + 1;
          right_exc_count[j] = right_exc_count[j] + 1;

          if (right_exc_count[i] == (PROBLEM_SIZE-3))
          {
            // Two edges must be included, and a node can't connect to itself, so 2 + 1 = 3
            right_include_edges(branches, right_inc_count, right_exc_count, i);
          }

          if (right_exc_count[j] == (PROBLEM_SIZE-3))
          {
            // Two edges must be included, and a node can't connect to itself, so 2 + 1 = 3
            right_include_edges(branches, right_inc_count, right_exc_count, j);
          }

          if(right_exc_count[i] > (PROBLEM_SIZE - 3) || right_exc_count[j] > (PROBLEM_SIZE - 3))
          {
            //printf("Error: right_exc_count > (NUM_NODES - 3)\n");
            delete branches;
            return nullptr;
          }

          done = true;
        }
      }
    }
  }

  return branches;
}

result_t
compute_lower_bound(int** graph, int explore[PROBLEM_SIZE][PROBLEM_SIZE])
{
  result_t sum = 0;

  for (size_t i = 0; i < PROBLEM_SIZE; ++i)
  {
    int min_edge1    = std::numeric_limits<int>::max();
    int min_edge2    = std::numeric_limits<int>::max();
    int forced_edge1 = std::numeric_limits<int>::max();
    int forced_edge2 = std::numeric_limits<int>::max();

    for(size_t j = 0; j < PROBLEM_SIZE; ++j)
    {
      if (i != j)
      {
        if(explore[i][j] == INCLUDE)
        {
          if (forced_edge1 == std::numeric_limits<int>::max())
          {
            forced_edge1 = graph[i][j];
            //printf("forced_edge1 = %d\n", graph[i][j]);
          }
          else if(forced_edge2 == std::numeric_limits<int>::max())
          {
            forced_edge2 = graph[i][j];
            //printf("forced_edge2 = %d\n", graph[i][j]);
          }
          else
          {
            //printf("Error: more than 2 edges selected\n");
          }
        }
        else if (explore[i][j] == UNDECIDED)
        {
          //A decision has not yet been made regarding this edge
          if (graph[i][j] < min_edge1)
          {
            min_edge2 = min_edge1;
            //printf("min_edge2 = %d\n", min_edge1);

            min_edge1 = graph[i][j];
            //printf("min_edge1 = %d\n", graph[i][j]);
          }
          else if(graph[i][j] < min_edge2)
          {
            min_edge2 = graph[i][j];
            //printf("min_edge2 = %d\n", graph[i][j]);
          }
        }
        //else {} //This edge will not be included in the explored path
      }
    }

    if (forced_edge1 != std::numeric_limits<int>::max() && forced_edge2 != std::numeric_limits<int>::max())
    {
      sum = sum + (result_t) (forced_edge1 + forced_edge2);
      //printf("Case 1: sum = %.2Lf\n", sum);
    }
    else if (forced_edge1 != std::numeric_limits<int>::max() && forced_edge2 == std::numeric_limits<int>::max())
    {
      sum = sum + (result_t) (forced_edge1 + min_edge1);
      //printf("Case 2: sum = %.2Lf\n", sum);
    }
    else if (forced_edge1 == std::numeric_limits<int>::max() && forced_edge2 == std::numeric_limits<int>::max())
    {
      sum = sum + (result_t) (min_edge1 + min_edge2);
      //printf("Case 3: sum = %.2Lf\n", sum);
    }
    else
    {
      //printf("Error: unexpected case in compute_lower_bound\n");
    }
  }

  return sum;
}

#if TBB_VERSION
template <class G, class T>
auto tsp_adaptive(G& taskgroup, T task) -> void // std::pair<D, size_t>
{
  while (true)
  {
    //base case
    if (task.left == true)
    {
      if (task.branch->left_branch[PROBLEM_SIZE-1][PROBLEM_SIZE-2] != UNDECIDED)
      {
        result_t result = compute_lower_bound(task.graph, task.branch->left_branch);
        //printf("task.res = %.2Lf\n", task.res);

        upd_lower_bound_if_needed(result);
        return;
      }
    }
    else
    {
      if (task.branch->right_branch[PROBLEM_SIZE-1][PROBLEM_SIZE-2] != UNDECIDED)
      {
        result_t result = compute_lower_bound(task.graph, task.branch->right_branch);

        upd_lower_bound_if_needed(result);
        return;
      }
    }

    BranchSet*           pbranch = task.left
                                   ? create_branch(task.branch->left_branch)
                                   : create_branch(task.branch->right_branch);

    if (pbranch == nullptr) return;

    tsp_task::branch_ptr branch(pbranch);

    T task2(task.graph, branch, true);
    taskgroup.run([&taskgroup, task2]()->void { tsp_adaptive(taskgroup, task2); });

    task = T{task.graph, branch, false};
  }
}

void tsp_launch(int** graph, BranchSet* branch, bool left)
{
  min_path.store(std::numeric_limits<result_t>::max());

  tbb::task_scheduler_init init(NUMTHREADS);
  tbb::task_group          g;
  tsp_task::branch_ptr    br(branch);

  tsp_adaptive(g, tsp_task(graph, br, left));
  g.wait();
  return;
}
#endif /* TBB_VERSION */

#if BLAZE_VERSION

struct Void
{
  Void operator+=(Void) { return *this; }
};

template <class T>
struct tsp_adaptive
{
  tsp_adaptive() {}

  Void operator()(uab::pool<T>& tasks, T task)
  {
    while(true)
    {
      //base case
      if(task.left == true)
      {
        if (task.branch->left_branch[PROBLEM_SIZE-1][PROBLEM_SIZE-2] != UNDECIDED)
        {
          result_t result = compute_lower_bound(task.graph, task.branch->left_branch);

          upd_lower_bound_if_needed(result);
          return Void();
        }
      }
      else
      {
        if (task.branch->right_branch[PROBLEM_SIZE-1][PROBLEM_SIZE-2] != UNDECIDED)
        {
          result_t result = compute_lower_bound(task.graph, task.branch->right_branch);

          upd_lower_bound_if_needed(result);
          return Void();
        }
      }

      BranchSet*           pbranch = task.left
                                     ? create_branch(task.branch->left_branch)
                                     : create_branch(task.branch->right_branch);

      if (pbranch == nullptr) return Void();

      tsp_task::branch_ptr branch(pbranch);

      tasks.enq(T{task.graph, branch, true});
      task = T{task.graph, branch, false};
    }
  }
};

void tsp_launch(int** graph, BranchSet* branch, bool left)
{
  min_path.store(std::numeric_limits<result_t>::max());

  tsp_adaptive<tsp_task>   fun;
  tsp_task::branch_ptr    br(branch);

  uab::execute_tasks(NUMTHREADS, fun, tsp_task(graph, br, left));
}
#endif /* BLAZE__VERSION */

#if OMP_VERSION

template <class T>
auto tsp_adaptive(T task) -> void // std::pair<D, size_t>
{
  while(true)
  {
    //base case
    if (task.left == true)
    {
      if(task.branch->left_branch[PROBLEM_SIZE-1][PROBLEM_SIZE-2] != UNDECIDED)
      {
        result_t result = compute_lower_bound(task.graph, task.branch->left_branch);

        upd_lower_bound_if_needed(result);
        return;
      }
    }
    else
    {
      if(task.branch->right_branch[PROBLEM_SIZE-1][PROBLEM_SIZE-2] != UNDECIDED)
      {
        result_t result = compute_lower_bound(task.graph, task.branch->right_branch);

        upd_lower_bound_if_needed(result);
        return;
      }
    }

    BranchSet*           pbranch = task.left
                                   ? create_branch(task.branch->left_branch)
                                   : create_branch(task.branch->right_branch);

    if (pbranch == nullptr) return;

    tsp_task::branch_ptr branch(pbranch);

    T task2(task.graph, branch, true);

    #pragma omp task firstprivate(task2)
    tsp_adaptive(task2);

    task = T{task.graph, branch, false};
  }
}

void tsp_launch(int** graph, BranchSet* branch, bool left)
{
  min_path.store(std::numeric_limits<result_t>::max());
  omp_set_num_threads(NUMTHREADS);

  tsp_task::branch_ptr    br(branch);

  #pragma omp parallel firstprivate(graph, br, left)
  #pragma omp single
  {
    #pragma omp taskgroup
    {
      tsp_adaptive(tsp_task(graph, br, left));
    }
  }
}

#endif /* OMP_VERSION */

#if CILK_VERSION

struct cilk_task
{
  int**      graph;
  BranchSet* branch;
  bool       left;

  cilk_task()
  : graph(nullptr), branch(nullptr), left(false)
  {}

  cilk_task(int** g, BranchSet* b, bool l)
  : graph(g), branch(b), left(l)
  {}
};


void set_cilk_workers(int n)
{
  assert(n <= 9999);

  char str[5];

  sprintf(str, "%d", n);

  bool success = __cilkrts_set_param("nworkers", str) != 0;
  assert(success);
}

void tsp_adaptive(cilk_task task);

static
BranchSet* tsp_adaptive_aux(cilk_task task)
{
  //base case
  if (task.left == true)
  {
    if(task.branch->left_branch[PROBLEM_SIZE-1][PROBLEM_SIZE-2] != UNDECIDED)
    {
      result_t result = compute_lower_bound(task.graph, task.branch->left_branch);

      upd_lower_bound_if_needed(result);
      return nullptr;
    }
  }
  else
  {
    if(task.branch->right_branch[PROBLEM_SIZE-1][PROBLEM_SIZE-2] != UNDECIDED)
    {
      result_t result = compute_lower_bound(task.graph, task.branch->right_branch);

      upd_lower_bound_if_needed(result);
      return nullptr;
    }
  }

  BranchSet*            branch = task.left
                                 ? create_branch(task.branch->left_branch)
                                 : create_branch(task.branch->right_branch);

  if (branch == nullptr) return nullptr;

  // Parallel Version
  cilk_task task2(task.graph, branch, true);

  cilk_spawn tsp_adaptive(task2);

  tsp_adaptive(cilk_task{task.graph, branch, false});

  return branch;
}

/// \brief
///   auxiliary function to deallocate BranchSet objects
/// \details
///   tsp_adaptive_aux creates a new BranchSet objects and passes it to its sub-task.
///   Since tsp_adaptive_aux implicitly waits until its sub-task has completed,
///   tsp_adaptive can free the allocated memory after its return.
void tsp_adaptive(cilk_task task)
{
  BranchSet* branch = tsp_adaptive_aux(task);

  if (branch) delete branch;
}


void tsp_launch(int** graph, BranchSet* branch, bool left)
{
  min_path.store(std::numeric_limits<result_t>::max());
  set_cilk_workers(NUMTHREADS);

  tsp_adaptive(cilk_task(graph, branch, left));

  if (branch) delete branch;
}
#endif /* CILK_VERSION */


#if QTHREADS_VERSION

void init_qthreads()
{
  std::stringstream str;

  str << "QTHREAD_HWPAR=" << NUMTHREADS;

  char* envset = new char[str.str().size()+1];

  memcpy(envset, str.str().c_str(), str.str().size()+1);
  putenv(envset);

  qthread_initialize();
}


struct qtask
{
  typedef tsp_task::branch_ptr branch_ptr;

  int**      graph;
  branch_ptr branch;
  bool       left;
  qt_sinc_t* sinc;

  qtask()
  : graph(nullptr), branch(nullptr), left(false), sinc(nullptr)
  {}

  qtask(int** g, branch_ptr b, bool l, qt_sinc_t* s)
  : graph(g), branch(b), left(l), sinc(s)
  {}
};


aligned_t tsp_adaptive(void* qtsk)
{
  qtask& task = *reinterpret_cast<qtask*>(qtsk);

  while(true)
  {
    //base case
    if (task.left == true)
    {
      if(task.branch->left_branch[PROBLEM_SIZE-1][PROBLEM_SIZE-2] != UNDECIDED)
      {
        result_t result = compute_lower_bound(task.graph, task.branch->left_branch);

        upd_lower_bound_if_needed(result);
        delete &task;
        qt_sinc_submit(task.sinc, nullptr);
        return aligned_t();
      }
    }
    else
    {
      if(task.branch->right_branch[PROBLEM_SIZE-1][PROBLEM_SIZE-2] != UNDECIDED)
      {
        result_t result = compute_lower_bound(task.graph, task.branch->right_branch);

        upd_lower_bound_if_needed(result);
        delete &task;
        qt_sinc_submit(task.sinc, nullptr);
        return aligned_t();
      }
    }

    BranchSet*           pbranch = task.left
                                   ? create_branch(task.branch->left_branch)
                                   : create_branch(task.branch->right_branch);

    if (pbranch == nullptr)
    {
      delete &task;
      qt_sinc_submit(task.sinc, nullptr);
      return aligned_t();
    }

    tsp_task::branch_ptr branch(pbranch);

    // Parallel Version
    qt_sinc_expect(task.sinc, 1);
    qthread_fork( tsp_adaptive,
                  new qtask(task.graph, branch, true, task.sinc),
                  nullptr
                );

    task = qtask{task.graph, branch, false, task.sinc};
  }
}


void tsp_launch(int** graph, BranchSet* branch, bool left)
{
  min_path.store(std::numeric_limits<result_t>::max());

  tsp_task::branch_ptr br(branch);
  qt_sinc_t*           sinc = qt_sinc_create(0, nullptr, nullptr, 0);

  tsp_adaptive(new qtask(graph, br, left, sinc));
  qt_sinc_wait(sinc, nullptr);
}

#endif /* QTHREADS_VERSION */


int main()
{
#if QTHREADS_VERSION
  init_qthreads();
#endif /* QTHREADS_VERSION */

  //input graph
  int graph[MAX_NODES][MAX_NODES] =
  {
    {0, 3, 4, 2, 7, 1, 1, 1, 1, 1, 1, 1},
    {3, 0, 4, 6, 3, 2, 2, 2, 2, 2, 2, 2},
    {4, 4, 0, 5, 8, 3, 3, 3, 3, 3, 3, 3},
    {2, 6, 5, 0, 6, 4, 4, 4, 4, 4, 4, 4},
    {7, 3, 8, 6, 0, 5, 5, 5, 5, 5, 5, 5},
    {1, 2, 3, 4, 5, 0, 6, 6, 6, 6, 6, 6},
    {1, 2, 3, 4, 5, 6, 0, 7, 7, 7, 7, 7},
    {1, 2, 3, 4, 5, 6, 7, 0, 8, 8, 8, 8},
    {1, 2, 3, 4, 5, 6, 7, 8, 0, 9, 9, 9},
    {1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 10, 10},
    {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 0, 11},
    {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0}
  };

  BranchSet* branch = new BranchSet;

  int explore[PROBLEM_SIZE][PROBLEM_SIZE];

  int** graph_ptr = (int**) malloc(PROBLEM_SIZE*sizeof(int*));

  for(size_t i = 0; i < PROBLEM_SIZE; ++i)
  {
    graph_ptr[i] = (int*) malloc(PROBLEM_SIZE*sizeof(int));
    for(size_t j = 0; j < PROBLEM_SIZE; ++j)
    {
      explore[i][j] = UNDECIDED;
      graph_ptr[i][j] = graph[i][j];
      branch->left_branch[i][j] = explore[i][j];
    }
  }

  typedef std::chrono::time_point<std::chrono::system_clock> time_point;
  time_point     starttime = std::chrono::system_clock::now();

  tsp_launch(graph_ptr, branch, true);

  time_point     endtime = std::chrono::system_clock::now();
  int            elapsedtime = std::chrono::duration_cast<std::chrono::milliseconds>(endtime-starttime).count();

  std::cout << "time = " << elapsedtime << "ms" << std::endl;
  std::cout << "Final result = " << min_path.load() << std::endl;
  std::cerr << elapsedtime << std::endl;

  for(size_t i = 0; i < PROBLEM_SIZE; ++i)
  {
    free(graph_ptr[i]);
  }

  free(graph_ptr);
  return 0;
}