/*
 * Code modified to run with the BLAZE task framework
 *
 * Peter Pirkelbauer, 2018
 * UAB - University of Alabama at Birmingham
 */


/**********************************************************************************************/
/*  This program is part of the Barcelona OpenMP Tasks Suite                                  */
/*  Copyright (C) 2009 Barcelona Supercomputing Center - Centro Nacional de Supercomputacion  */
/*  Copyright (C) 2009 Universitat Politecnica de Catalunya                                   */
/*                                                                                            */
/*  This program is free software; you can redistribute it and/or modify                      */
/*  it under the terms of the GNU General Public License as published by                      */
/*  the Free Software Foundation; either version 2 of the License, or                         */
/*  (at your option) any later version.                                                       */
/*                                                                                            */
/*  This program is distributed in the hope that it will be useful,                           */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of                            */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                             */
/*  GNU General Public License for more details.                                              */
/*                                                                                            */
/*  You should have received a copy of the GNU General Public License                         */
/*  along with this program; if not, write to the Free Software                               */
/*  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA            */
/**********************************************************************************************/

/*
 * Original code from the Cilk project
 *
 * Copyright (c) 2000 Massachusetts Institute of Technology
 * Copyright (c) 2000 Matteo Frigo
 */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>

#include <atomic>
#include <chrono>

#include "../common/bots.hpp"

#define MAX_ITEMS 256

#ifndef NUMTHREADS
#define NUMTHREADS (20)
#endif /* NUMTHREADS */

#if OMP_VERSION
#include <omp.h>
#endif /* OMP_VERSION */

#if BOTS_VERSION
// close to the original version in the Barcelona OpenMP Testing Suite (BOTS).
#include <omp.h>
#endif /* BOTS_VERSION */

#if TBB_VERSION
#include <mutex>
#include <tbb/task_scheduler_init.h>
#include <tbb/task_group.h>
#endif /* TBB_VERSION */

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
#include <qthread/qthread.hpp>
#include <qthread/sinc.h>
#include "../common/qthreads.hpp"
#endif /* QTHREADS_VERSION */


struct item
{
  int value;
  int weight;
};

std::atomic<int> best_so_far; //< new best so far

int bots_best_so_far;         //< original best so far -- only for BOTS code

int compare(struct item *a, struct item *b)
{
     double c = ((double) a->value / a->weight) -
     ((double) b->value / b->weight);

     if (c > 0) return -1;
     if (c < 0) return 1;
     return 0;
}

int read_input(const char *filename, struct item *items, int *capacity, int *n)
{
     int i;
     FILE *f;

     if (filename == NULL) filename = "\0";
     f = fopen(filename, "r");
     if (f == NULL)
     {
       fprintf(stderr, "open_input(\"%s\") failed\n", filename);
       exit(1);
     }
     /* format of the input: #items capacity value1 weight1 ... */
     fscanf(f, "%d", n);
     fscanf(f, "%d", capacity);

     for (i = 0; i < *n; ++i)
       fscanf(f, "%d %d", &items[i].value, &items[i].weight);

     fclose(f);

     /* sort the items on decreasing order of value/weight */
     /* cilk2c is fascist in dealing with pointers, whence the ugly cast */
     qsort(items, *n, sizeof(struct item), (int (*)(const void *, const void *)) compare);
     return 0;
}

void update_best_so_far(int val)
{
  int curr = best_so_far.load(std::memory_order_relaxed);

  while (val > curr)
  {
    best_so_far.compare_exchange_weak(curr, val, std::memory_order_relaxed, std::memory_order_relaxed);
  }
}

struct knapsack_task
{
  item* e;
  int   c;
  int   n;
  int   v;
};


#if BLAZE_VERSION

struct knapsack_par
{
  template <class P>
  auto operator()(P& pool, knapsack_task task) -> uab::Void
  {
    for (;;)
    {
      /* base case: full knapsack or no items */
      if (task.c < 0) return uab::Void();

      /* feasible solution, with value v */
      if (task.n == 0 || task.c == 0)
      {
        update_best_so_far(task.v);
        return uab::Void();
      }

      double ub = (double) task.v + task.c * task.e->value / task.e->weight;

      // prune if it is worse than the best available alternative
      if (ub < best_so_far.load(std::memory_order_relaxed)) return uab::Void();

      /* compute the best solution without the current item in the knapsack */
      // #pragma omp task untied firstprivate(e,c,n,v,l) shared(without)
      pool.enq(knapsack_task{ task.e + 1, task.c, task.n - 1, task.v });

      /* compute the best solution with the current item in the knapsack */
      // #pragma omp task untied firstprivate(e,c,n,v,l) shared(with)
      task = knapsack_task{ task.e + 1, task.c - task.e->weight, task.n - 1, task.v + task.e->value };
    }
  }
};

auto knapsack(item *e, int c, int n) -> int
{
  uab::execute_tasks(NUMTHREADS, knapsack_par(), knapsack_task{e, c, n, 0} );

  return best_so_far.load(std::memory_order_relaxed);
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


auto knapsack_par(knapsack_task task) -> void
{
  for (;;)
  {
    /* base case: full knapsack or no items */
    if (task.c < 0) return;

    /* feasible solution, with value v */
    if (task.n == 0 || task.c == 0)
    {
      update_best_so_far(task.v);
      return;
    }

    double ub = (double) task.v + task.c * task.e->value / task.e->weight;

    // prune if it is worse than the best available alternative
    if (ub < best_so_far.load(std::memory_order_relaxed)) return;

    /* compute the best solution without the current item in the knapsack */
    // #pragma omp task untied firstprivate(e,c,n,v,l) shared(without)
    cilk_spawn knapsack_par(knapsack_task{ task.e + 1, task.c, task.n - 1, task.v });

    /* compute the best solution with the current item in the knapsack */
    // #pragma omp task untied firstprivate(e,c,n,v,l) shared(with)
    task = knapsack_task{ task.e + 1, task.c - task.e->weight, task.n - 1, task.v + task.e->value };
  };
}

auto knapsack(item *e, int c, int n) -> int
{
  set_cilk_workers(NUMTHREADS);
  knapsack_par(knapsack_task{e, c, n, 0});

  return best_so_far.load(std::memory_order_relaxed);
}

#endif /* CILK_VERSION */

#if OMP_VERSION

auto knapsack_par(knapsack_task task) -> void
{
  for (;;)
  {
    /* base case: full knapsack or no items */
    if (task.c < 0) return;

    /* feasible solution, with value v */
    if (task.n == 0 || task.c == 0)
    {
      update_best_so_far(task.v);
      return;
    }

    double ub = (double) task.v + task.c * task.e->value / task.e->weight;

    // prune if it is worse than the best available alternative
    if (ub < best_so_far.load(std::memory_order_relaxed)) return;

    /* compute the best solution without the current item in the knapsack */
    #pragma omp task untied firstprivate(task)
    knapsack_par(knapsack_task{ task.e + 1, task.c, task.n - 1, task.v });

    /* compute the best solution with the current item in the knapsack */
    task = knapsack_task{ task.e + 1, task.c - task.e->weight, task.n - 1, task.v + task.e->value };
  }
}

auto knapsack (item* e, int c, int n) -> int
{
     omp_set_num_threads(NUMTHREADS);

     best_so_far = INT_MIN;

     #pragma omp parallel firstprivate(e,c,n)
     {
        #pragma omp single
        #pragma omp taskgroup
        {
           knapsack_par(knapsack_task{e, c, n, 0});
        }
     }

     return best_so_far.load(std::memory_order_relaxed);
}

#endif /* OMP_VERSION */

#if BOTS_VERSION

/*
 * return the optimal solution for n items (first is e) and
 * capacity c. Value so far is v.
 */
void knapsack_par(struct item *e, int c, int n, int v, int *sol, int l)
{
     int with, without, best;
     double ub;

     /* base case: full knapsack or no items */
     if (c < 0)
     {
         *sol = INT_MIN;
         return;
     }

     /* feasible solution, with value v */
     if (n == 0 || c == 0)
     {
         *sol = v;
         return;
     }

     ub = (double) v + c * e->value / e->weight;

     if (ub < bots_best_so_far) { /* prune ! */
          *sol = INT_MIN;
          return;
     }

     /* compute the best solution without the current item in the knapsack */
     #pragma omp task untied firstprivate(e,c,n,v,l) shared(without)
     knapsack_par(e + 1, c, n - 1, v, &without, l+1);

     /* compute the best solution with the current item in the knapsack */
     #pragma omp task untied firstprivate(e,c,n,v,l) shared(with)
     knapsack_par(e + 1, c - e->weight, n - 1, v + e->value, &with, l+1);

     #pragma omp taskwait
     best = with > without ? with : without;

     /*
      * notice the race condition here. The program is still
      * correct, in the sense that the best solution so far
      * is at least bots_best_so_far. Moreover bots_best_so_far gets updated
      * when returning, so eventually it should get the right
      * value. The program is highly non-deterministic.
      */
     if (best > bots_best_so_far) bots_best_so_far = best;

     *sol = best;
}

int knapsack (struct item *e, int c, int n)
{
     int sol = INT_MIN;

     omp_set_num_threads(NUMTHREADS);
     bots_best_so_far = INT_MIN;

     #pragma omp parallel
     {
        #pragma omp single
        #pragma omp task untied
        {
           knapsack_par(e, c, n, 0, &sol, 0);
        }
     }

     return sol;
}

#endif /* BOTS_VERSION */


#if TBB_VERSION

template <class G>
auto knapsack_par(G& taskgroup, knapsack_task task) -> void
{
  for (;;)
  {
    /* base case: full knapsack or no items */
    if (task.c < 0) return;

    /* feasible solution, with value v */
    if (task.n == 0 || task.c == 0)
    {
      update_best_so_far(task.v);
      return;
    }

    double ub = (double) task.v + task.c * task.e->value / task.e->weight;

    // prune if it is worse than the best available alternative
    if (ub < best_so_far.load(std::memory_order_relaxed)) return;

    /* compute the best solution without the current item in the knapsack */
    // #pragma omp task untied firstprivate(e,c,n,v,l) shared(without)
    taskgroup.run( [&taskgroup, task]() -> void
                   {
                     knapsack_par(taskgroup, knapsack_task{ task.e + 1, task.c, task.n - 1, task.v });
                   }
                 );

    /* compute the best solution with the current item in the knapsack */
    // #pragma omp task untied firstprivate(e,c,n,v,l) shared(with)
    task = knapsack_task{ task.e + 1, task.c - task.e->weight, task.n - 1, task.v + task.e->value };
  }
}

auto knapsack(item *e, int c, int n) -> int
{
   tbb::task_scheduler_init init(NUMTHREADS);
   tbb::task_group          g;

   knapsack_par(g, knapsack_task{e, c, n, 0});
   g.wait();

  return best_so_far.load(std::memory_order_relaxed);
}


#endif /* TBB_VERSION */


#if QTHREADS_VERSION

struct qtask
{
  item*      e;
  int        c;
  int        n;
  int        v;
  qt_sinc_t* sinc;
};


aligned_t knapsack_par(void* qtsk)
{
  qtask& task = * reinterpret_cast<qtask*>(qtsk);

  for (;;)
  {
    /* base case: full knapsack or no items */
    if (task.c < 0)
    {
      qt_sinc_submit(task.sinc, nullptr);
      delete &task;
      return aligned_t();
    }

    /* feasible solution, with value v */
    if (task.n == 0 || task.c == 0)
    {
      update_best_so_far(task.v);
      qt_sinc_submit(task.sinc, nullptr);
      delete &task;
      return aligned_t();
    }

    double ub = (double) task.v + task.c * task.e->value / task.e->weight;

    // prune if it is worse than the best available alternative
    if (ub < best_so_far.load(std::memory_order_relaxed))
    {
      qt_sinc_submit(task.sinc, nullptr);
      delete &task;
      return aligned_t();
    }

    /* compute the best solution without the current item in the knapsack */
    qt_sinc_expect(task.sinc, 1);
    qthread_fork( knapsack_par,
                  new qtask{ task.e + 1, task.c, task.n - 1, task.v, task.sinc },
                  nullptr
                );

    /* compute the best solution with the current item in the knapsack */
    task = qtask{ task.e + 1,
                  task.c - task.e->weight,
                  task.n - 1,
                  task.v + task.e->value,
                  task.sinc
                };
  }
}

auto knapsack(item *e, int c, int n) -> int
{
  qt_sinc_t* sinc   = qt_sinc_create(0, nullptr, nullptr, 1);

  knapsack_par(new qtask{e, c, n, 0, sinc});

  qt_sinc_wait(sinc, nullptr);

  return best_so_far.load(std::memory_order_relaxed);
}

#endif /* QTHREADS_VERSION */




void knapsack_seq(struct item *e, int c, int n, int v, int *sol)
{
     int with, without, best;
     double ub;

     /* base case: full knapsack or no items */
     if (c < 0)
     {
         *sol = INT_MIN;
         return;
     }

     /* feasible solution, with value v */
     if (n == 0 || c == 0)
     {
         *sol = v;
         return;
     }

     ub = (double) v + c * e->value / e->weight;

     if (ub < bots_best_so_far) {
      /* prune ! */
          *sol = INT_MIN;
          return;
     }
     /*
      * compute the best solution without the current item in the knapsack
      */
     knapsack_seq(e + 1, c, n - 1, v, &without);

     /* compute the best solution with the current item in the knapsack */
     knapsack_seq(e + 1, c - e->weight, n - 1, v + e->value, &with);

     best = with > without ? with : without;

     /*
      * notice the race condition here. The program is still
      * correct, in the sense that the best solution so far
      * is at least best_so_far. Moreover best_so_far gets updated
      * when returning, so eventually it should get the right
      * value. The program is highly non-deterministic.
      */
     if (best > bots_best_so_far) bots_best_so_far = best;

     *sol = best;
}

int knapsack_main_seq(struct item *e, int c, int n)
{
  int sol = INT_MIN;

  bots_best_so_far = INT_MIN;

  knapsack_seq(e, c, n, 0, &sol);
  return sol;
}

int knapsack_check (int sol1, int sol2)
{
   if (sol1 == sol2) return BOTS_RESULT_SUCCESSFUL;

   return BOTS_RESULT_UNSUCCESSFUL;
}

int main(int argc, char** argv)
{
  typedef std::chrono::time_point<std::chrono::system_clock> time_point;

  std::string inputfile = "data/knapsack-044.input";
  item        items[MAX_ITEMS];
  int         n;
  int         capacity;

  if (argc > 1) inputfile = argv[1];

  std::cout << "loading " << inputfile << std::endl;
  read_input(inputfile.c_str(), items, &capacity, &n);

#if QTHREADS_VERSION
  init_qthreads(NUMTHREADS);
#endif /* QTHREADS_VERSION */

  time_point starttime = std::chrono::system_clock::now();
  int        solpar = knapsack(items, capacity, n);
  time_point endtime = std::chrono::system_clock::now();
  int        elapsedtime = std::chrono::duration_cast<std::chrono::milliseconds>(endtime-starttime).count();
/*
  int        solseq = knapsack_main_seq(items, capacity, n);
  int        valid = knapsack_check(solseq, solpar);

  if (valid 1= BOTS_RESULT_SUCCESSFUL)
  {
    // print timing
    std::cout << "result mismatch" << std::endl;
    std::cerr << "result mismatch" << std::endl;
    return 1;
  }
*/

  // print timing
  std::cout << "time = " << elapsedtime << "ms; solution = " << solpar << std::endl;
  std::cerr << elapsedtime << std::endl;

  return 0;
}
