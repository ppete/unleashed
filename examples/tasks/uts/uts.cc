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
/**********************************************************************************************/
/*
 * Copyright (c) 2007 The Unbalanced Tree Search (UTS) Project Team:
 * -----------------------------------------------------------------
 *
 *  This file is part of the unbalanced tree search benchmark.  This
 *  project is licensed under the MIT Open Source license.  See the LICENSE
 *  file for copyright and licensing information.
 *
 *  UTS is a collaborative project between researchers at the University of
 *  Maryland, the University of North Carolina at Chapel Hill, and the Ohio
 *  State University.
 *
 * University of Maryland:
 *   Chau-Wen Tseng(1)  <tseng at cs.umd.edu>
 *
 * University of North Carolina, Chapel Hill:
 *   Jun Huan         <huan,
 *   Jinze Liu         liu,
 *   Stephen Olivier   olivier,
 *   Jan Prins*        prins at cs.umd.edu>
 *
 * The Ohio State University:
 *   James Dinan      <dinan,
 *   Gerald Sabin      sabin,
 *   P. Sadayappan*    saday at cse.ohio-state.edu>
 *
 * Supercomputing Research Center
 *   D. Pryor
 *
 * (1) - indicates project PI
 *
 * UTS Recursive Depth-First Search (DFS) version developed by James Dinan
 *
 * Adapted for OpenMP 3.0 Task-based version by Stephen Olivier
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 *
 */
#include <iostream>
#include <iomanip>
#include <chrono>
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <list>
#include <vector>

#include "../common/common-includes.hpp"
#include "../common/bots.hpp"
#include "atomicutil.hpp"

//~ #include "app-desc.h"
//~ #include "bots.h"
#include "uts.h"

typedef unsigned long long count_t;

/***********************************************************
 *  Global state                                           *
 ***********************************************************/
count_t  nLeaves = 0;
int      maxTreeDepth = 0;

/***********************************************************
 * Tree generation strategy is controlled via various      *
 * parameters set from the command line.  The parameters   *
 * and their default values are given below.               *
 * Trees are generated using a Galton-Watson process, in   *
 * which the branching factor of each node is a random     *
 * variable.                                               *
 *                                                         *
 * The random variable follow a binomial distribution.     *
 ***********************************************************/
double   b_0   = 4.0; // default branching factor at the root
int      rootId = 0;   // default seed for RNG state at root

/***********************************************************
 *  The branching factor at the root is specified by b_0.
 *  The branching factor below the root follows an
 *     identical binomial distribution at all nodes.
 *  A node has m children with prob q, or no children with
 *     prob (1-q).  The expected branching factor is q * m.
 *
 *  Default parameter values
 ***********************************************************/
int      nonLeafBF   = 4;            // m
double   nonLeafProb = 15.0 / 64.0;  // q

/***********************************************************
 * compute granularity - number of rng evaluations per
 * tree node
 ***********************************************************/
int      computeGranularity = 1;

/***********************************************************
 * expected results for execution
 ***********************************************************/
count_t  exp_tree_size = 0;
int      exp_tree_depth = 0;
count_t  exp_num_leaves = 0;

/***********************************************************
 *  FUNCTIONS                                              *
 ***********************************************************/

// Interpret 32 bit positive integer as value on [0,1)
double rng_toProb(int n)
{
  if (n < 0) {
    printf("*** toProb: rand n = %d out of range\n",n);
  }
  return ((n<0)? 0.0 : ((double) n)/2147483648.0);
}

void uts_initRoot(Node * root)
{
   root->height = 0;
   root->numChildren = -1;      // means not yet determined
   rng_init(root->state.state, rootId);

   bots_message("Root node at %p\n", root);
}


int uts_numChildren_bin(Node * parent)
{
  // distribution is identical everywhere below root
  int    v = rng_rand(parent->state.state);
  double d = rng_toProb(v);

  return (d < nonLeafProb) ? nonLeafBF : 0;
}

int uts_numChildren(Node *parent)
{
  int numChildren = 0;

  /* Determine the number of children */
  if (parent->height == 0) numChildren = (int) floor(b_0);
  else numChildren = uts_numChildren_bin(parent);

  // limit number of children
  // only a BIN root can have more than MAXNUMCHILDREN
  if (parent->height == 0) {
    int rootBF = (int) ceil(b_0);
    if (numChildren > rootBF) {
      bots_debug("*** Number of children of root truncated from %d to %d\n", numChildren, rootBF);
      numChildren = rootBF;
    }
  }
  else {
    if (numChildren > MAXNUMCHILDREN) {
      bots_debug("*** Number of children truncated from %d to %d\n", numChildren, MAXNUMCHILDREN);
      numChildren = MAXNUMCHILDREN;
    }
  }

  return numChildren;
}

/***********************************************************
 * Recursive depth-first implementation                    *
 ***********************************************************/

#if BOTS_VERSION

count_t parTreeSearch(int depth, Node *parent, int numChildren)
{
  Node             children[numChildren];
  count_t          partialCount[numChildren];

  // Recurse on the children
  for (int i = 0; i < numChildren; i++) {
     Node* child = &children[i];

     child->height = parent->height + 1;

     // The following line is the work (one or more SHA-1 ops)
     for (int j = 0; j < computeGranularity; j++) {
        rng_spawn(parent->state.state, child->state.state, i);
     }

     child->numChildren = uts_numChildren(child);

     #pragma omp task untied firstprivate(depth, i, child) shared(partialCount)
     partialCount[i] = parTreeSearch(depth+1, child, child->numChildren);
  }

  #pragma omp taskwait

  count_t subtreesize = 1;

  for (int i = 0; i < numChildren; i++) {
     subtreesize += partialCount[i];
  }

  return subtreesize;
}


count_t parallel_uts(Node* root, size_t numthreads)
{
   omp_set_num_threads(numthreads);

   count_t num_nodes = 0 ;
   root->numChildren = uts_numChildren(root);

   #pragma omp parallel
      #pragma omp single nowait
      #pragma omp task untied
      num_nodes = parTreeSearch( 0, root, root->numChildren );

   return num_nodes;
}
#endif /* BOTS_VERSION */

#if OMP_VERSION

typedef uab::aligned_type<count_t, CACHELINESZ> accumulator_t;

static inline
size_t thread_num()
{
  return omp_get_thread_num();
}

void parTreeSearch(int depth, Node* parent, int numChildren, accumulator_t* sum)
{
  // Recurse on the children
  for (int i = 0; i < numChildren; ++i)
  {
     Node* nodePtr = new Node;

     nodePtr->height = parent->height + 1;

     // The following line is the work (one or more SHA-1 ops)
     for (int j = 0; j < computeGranularity; ++j)
     {
       rng_spawn(parent->state.state, nodePtr->state.state, i);
     }

     nodePtr->numChildren = uts_numChildren(nodePtr);

     #pragma omp task untied firstprivate(depth, nodePtr, sum)
     parTreeSearch(depth+1, nodePtr, nodePtr->numChildren, sum);
  }

  delete parent;
  ++sum[thread_num()].val;
}


count_t parallel_uts (Node* root, size_t numthreads)
{
   omp_set_num_threads(NUMTHREADS);

   count_t num_nodes = 0;
   accumulator_t      sum[NUMTHREADS];

   root->numChildren = uts_numChildren(root);

   #pragma omp parallel firstprivate(root) shared(sum)
   {
      #pragma omp single nowait
      #pragma omp taskgroup
      parTreeSearch( 0, root, root->numChildren, sum );
   }

   while (numthreads) num_nodes += sum[--numthreads].val;

   return num_nodes;
}
#endif /* OMP_VERSION */

#if BLAZE_VERSION

struct par_tree_task
{
  int   depth;
  Node* parent;
  int   numChildren;
};

struct par_tree_search
{
  template <class P>
  count_t
  operator()(P& pool, par_tree_task task)
  {
    // Recurse on the children
    for (int i = 0; i < task.numChildren; ++i)
    {
       Node* nodePtr = new Node;

       nodePtr->height = task.parent->height + 1;

       // The following line is the work (one or more SHA-1 ops)
       for (int j = 0; j < computeGranularity; ++j)
       {
         rng_spawn(task.parent->state.state, nodePtr->state.state, i);
       }

       nodePtr->numChildren = uts_numChildren(nodePtr);

       pool.enq(par_tree_task{task.depth+1, nodePtr, nodePtr->numChildren});
    }

    delete task.parent;
    return 1;
  }
};


count_t parallel_uts (Node* root, size_t numthreads)
{
   root->numChildren = uts_numChildren(root);

   return uab::execute_tasks( numthreads,
                              par_tree_search(),
                              par_tree_task { 0, root, root->numChildren }
                            );
}
#endif /* BLAZE_VERSION */


#if CILK_VERSION

void cilk_init(int workers, std::string stacksz)
{
  assert(workers <= 9999);

  char str[5];

  sprintf(str, "%d", workers);

  bool stackset = __cilkrts_set_param("stack size", stacksz.c_str()) == 0;
  assert(stackset);

  bool success = __cilkrts_set_param("nworkers", str) == 0;
  assert(success);
}

void
cilk_tree_search(int depth, Node* parent, int numChildren, cilk::reducer_opadd<count_t>& sum, Node* children);

void
par_tree_search(int depth, Node* parent, int numChildren, cilk::reducer_opadd<count_t>& sum)
{
  Node children[numChildren];

  cilk_tree_search(depth, parent, numChildren, sum, children);
}

void
cilk_tree_search(int depth, Node* parent, int numChildren, cilk::reducer_opadd<count_t>& sum, Node* children)
{
  // Recurse on the children
  for (int i = 0; i < numChildren; ++i)
  {
    Node* child = &children[i];

    child->height = parent->height + 1;

    // The following line is the work (one or more SHA-1 ops)
    for (int j = 0; j < computeGranularity; ++j)
    {
      rng_spawn(parent->state.state, child->state.state, i);
    }

    child->numChildren = uts_numChildren(child);
    cilk_spawn par_tree_search(depth+1, child, child->numChildren, sum);
  }

  sum += 1;
}


count_t parallel_uts (Node* root, size_t numthreads)
{
  cilk_init(numthreads, "4000000");

  cilk::reducer_opadd<count_t> sum(0);

  root->numChildren = uts_numChildren(root);
  par_tree_search(0, root, root->numChildren, sum);

  return sum.get_value();
}
#endif /* CILK_VERSION */


#if QTHREADS_VERSION

// \note code adopted from the QTHREADS benchmark

/******************************************************
* Unbalanced Tree Search v2.1                         *
* Based on the implementation available at            *
*     http://sourceforge.net/projects/uts-benchmark   *
******************************************************/

static qt_sinc_t* sinc = nullptr;

// Notes:
// -    Each task receives distinct copy of parent
// -    Copy of child is shallow, be careful with `state` member
static aligned_t visit(void *args_)
{
  Node*  parent        = (Node*)args_;
  int    parent_height = parent->height;
  int    num_children  = parent->numChildren;
  node_t child;

  qt_sinc_expect(sinc, num_children);

  // Spawn children, if any
  for (int i = 0; i < num_children; ++i)
  {
    child.height = parent_height + 1;

    // The following line is the work (one or more SHA-1 ops)
    for (int j = 0; j < computeGranularity; ++j)
    {
      rng_spawn(parent->state.state, child.state.state, i);
    }

    child.numChildren = uts_numChildren(&child);
    qthread_fork_syncvar_copyargs(visit, &child, sizeof(node_t), NULL);
  }

  count_t value = 1;
  qt_sinc_submit(sinc, &value);
  return 0;
}

static
void my_incr(void *tgt, const void *src)
{
  *(count_t *)tgt += *(count_t *)src;
}


count_t parallel_uts (Node* root, size_t /* numthreads */)
{
  root->numChildren = uts_numChildren(root);

  count_t initial_value   = 0;
  count_t total_num_nodes = 0;

  sinc = qt_sinc_create(sizeof(count_t), &initial_value, my_incr, 1);
  visit(root);

  qt_sinc_wait(sinc, &total_num_nodes);
  qt_sinc_destroy(sinc);

  return total_num_nodes;
}


#endif /* QTHREADS_VERSION */

#if TBB_VERSION

template <class G>
void
par_tree_search(G& g, int depth, Node* parent, int numChildren, uab::simple_reducer<count_t>& sum)
{
  // Recurse on the children
  for (int i = 0; i < numChildren; ++i)
  {
    Node* child = new Node;

    child->height = parent->height + 1;

    // The following line is the work (one or more SHA-1 ops)
    for (int j = 0; j < computeGranularity; ++j)
    {
      rng_spawn(parent->state.state, child->state.state, i);
    }

    child->numChildren = uts_numChildren(child);
    g.run( [&g, depth, child, &sum]() -> void
           {
             par_tree_search(g, depth+1, child, child->numChildren, sum);
           }
         );
  }

  sum += 1;
  delete parent;
}


count_t parallel_uts(Node* root, size_t numthreads)
{
  tbb::task_scheduler_init     init(numthreads);
  tbb::task_group              g;
  uab::simple_reducer<count_t> reducer;

  root->numChildren = uts_numChildren(root);

  par_tree_search(g, 0, root, root->numChildren, reducer);

  g.wait();
  return reducer.get_value();
}

#endif /* TBB_VERSION */



void uts_read_file ( const char *filename )
{
   FILE *fin;

   if ((fin = fopen(filename, "r")) == NULL) {
      bots_message("Could not open input file (%s)\n", filename);
      exit (-1);
   }

   fscanf(fin,"%lf %lf %d %d %d %llu %d %llu",
             &b_0,
             &nonLeafProb,
             &nonLeafBF,
             &rootId,
             &computeGranularity,
             &exp_tree_size,
             &exp_tree_depth,
             &exp_num_leaves
   );
   fclose(fin);

   computeGranularity = max(1,computeGranularity);

   // Printing input data
   bots_message("\n");
   bots_message("Root branching factor                = %f\n", b_0);
   bots_message("Root seed (0 <= 2^31)                = %d\n", rootId);
   bots_message("Probability of non-leaf node         = %f\n", nonLeafProb);
   bots_message("Number of children for non-leaf node = %d\n", nonLeafBF);
   bots_message("E(n)                                 = %f\n", (double) ( nonLeafProb * nonLeafBF ) );
   bots_message("E(s)                                 = %f\n", (double) ( 1.0 / (1.0 - nonLeafProb * nonLeafBF) ) );
   bots_message("Compute granularity                  = %d\n", computeGranularity);
   bots_message("Random number generator              = "); rng_showtype();
}

#if 0 /* NOT_USED */
void uts_show_stats( void )
{
   int nPes = atoi(bots_resources);
   int chunkSize = 0;

   bots_message("\n");
   bots_message("Tree size                            = %llu\n", (count_t)  bots_number_of_tasks );
   bots_message("Maximum tree depth                   = %d\n", maxTreeDepth );
   bots_message("Chunk size                           = %d\n", chunkSize );
   bots_message("Number of leaves                     = %llu (%.2f%%)\n", nLeaves, nLeaves/(float)bots_number_of_tasks*100.0 );
   bots_message("Number of PE's                       = %.4d threads\n", nPes );
   bots_message("Wallclock time                       = %.3f sec\n", bots_time_program );
   bots_message("Overall performance                  = %.0f nodes/sec\n", (bots_number_of_tasks / bots_time_program) );
   bots_message("Performance per PE                   = %.0f nodes/sec\n", (bots_number_of_tasks / bots_time_program / nPes) );
}
#endif /* NOT_USED */

int uts_check_result ( count_t num_tasks )
{
  int answer = BOTS_RESULT_SUCCESSFUL;

  if ( num_tasks != exp_tree_size )
  {
    answer = BOTS_RESULT_UNSUCCESSFUL;
    bots_message("Incorrect tree size (%llu instead of %llu).\n", num_tasks, exp_tree_size);
  }

  return answer;
}


int main(int argc, char** argv)
{
  typedef std::chrono::time_point<std::chrono::system_clock> time_point;

  std::string inputfile = "data/medium.input";
  Node*       root = new Node;

  if (argc > 1) inputfile = argv[1];

  std::cout << "loading " << inputfile << std::endl;
  uts_read_file(inputfile.c_str());
  uts_initRoot(root);

#if QTHREADS_VERSION
  init_qthreads(NUMTHREADS);
#endif /* QTHREADS_VERSION */

  time_point starttime   = std::chrono::system_clock::now();
  count_t    solpar      = parallel_uts(root, NUMTHREADS);
  time_point endtime     = std::chrono::system_clock::now();
  int        elapsedtime = std::chrono::duration_cast<std::chrono::milliseconds>(endtime-starttime).count();

  if (uts_check_result(solpar) != BOTS_RESULT_SUCCESSFUL)
  {
    std::cout << "result mismatch" << std::endl;
    std::cerr << "result mismatch" << std::endl;
    return 1;
  }

  // print timing
  std::cout << "time = " << elapsedtime << "ms; solution = " << solpar << std::endl;
  std::cerr << elapsedtime << std::endl;

  return 0;
}
