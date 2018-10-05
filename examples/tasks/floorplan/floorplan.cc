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

/* Original code from the Application Kernel Matrix by Cray */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <mutex>
#include <atomic>
#include <memory>

#include "../common/bots.hpp"

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

#define ROWS 64
#define COLS 64
#define DMAX 64
#define max(a, b) ((a > b) ? a : b)
#define min(a, b) ((a < b) ? a : b)

int solution = -1;

typedef int  coor[2];
typedef char ibrd[ROWS][COLS];
typedef char (*pibrd)[COLS];

FILE * inputFile;

struct cell {
  int   n;
  coor *alt;
  int   top;
  int   bot;
  int   lhs;
  int   rhs;
  int   left;
  int   above;
  int   next;
};

typedef std::shared_ptr<cell> cell_ptr;
//~ typedef uab::counted_array<cell> cell_ptr;

static cell*             gcells;
static std::atomic<int>  MIN_AREA;
static ibrd              BEST_BOARD;
static coor              MIN_FOOTPRINT;
static int               N;

/* compute all possible locations for nw corner for cell */
static int starts(int id, int shape, coor *NWS, struct cell *cells) {
  int i, n, top, bot, lhs, rhs;
  int rows, cols, left, above;

/* size of cell */
  rows  = cells[id].alt[shape][0];
  cols  = cells[id].alt[shape][1];

/* the cells to the left and above */
  left  = cells[id].left;
  above = cells[id].above;

/* if there is a vertical and horizontal dependence */
  if ((left >= 0) && (above >= 0)) {

     top = cells[above].bot + 1;
     lhs = cells[left].rhs + 1;
     bot = top + rows;
     rhs = lhs + cols;

/* if footprint of cell touches the cells to the left and above */
     if ((top <= cells[left].bot) && (bot >= cells[left].top) &&
         (lhs <= cells[above].rhs) && (rhs >= cells[above].lhs))
          { n = 1; NWS[0][0] = top; NWS[0][1] = lhs;  }
     else { n = 0; }

/* if there is only a horizontal dependence */
   } else if (left >= 0) {

/* highest initial row is top of cell to the left - rows */
     top = max(cells[left].top - rows + 1, 0);
/* lowest initial row is bottom of cell to the left */
     bot = min(cells[left].bot, ROWS);
     n   = bot - top + 1;

     for (i = 0; i < n; i++) {
         NWS[i][0] = i + top;
         NWS[i][1] = cells[left].rhs + 1;
     }

  } else {

/* leftmost initial col is lhs of cell above - cols */
     lhs = max(cells[above].lhs - cols + 1, 0);
/* rightmost initial col is rhs of cell above */
     rhs = min(cells[above].rhs, COLS);
     n   = rhs - lhs + 1;

     for (i = 0; i < n; i++) {
         NWS[i][0] = cells[above].bot + 1;
         NWS[i][1] = i + lhs;
  }  }

  return (n);
}



/* lay the cell down on the board in the rectangular space defined
   by the cells top, bottom, left, and right edges. If the cell can
   not be layed down, return 0; else 1.
*/
static int lay_down(int id, ibrd board, struct cell *cells) {
  int  i, j, top, bot, lhs, rhs;

  top = cells[id].top;
  bot = cells[id].bot;
  lhs = cells[id].lhs;
  rhs = cells[id].rhs;

  for (i = top; i <= bot; i++) {
  for (j = lhs; j <= rhs; j++) {
      if (board[i][j] == 0) board[i][j] = (char)id;
      else                  return(0);
  } }

  return (1);
}


#define read_integer(file,var) \
  if ( fscanf(file, "%d", &var) == EOF ) {\
	bots_message(" Bogus input file\n");\
	exit(-1);\
  }

static void read_inputs() {
  int i, j, n;

  read_integer(inputFile,n);
  N = n;

  gcells = (struct cell *) malloc((n + 1) * sizeof(struct cell));

  gcells[0].n     =  0;
  gcells[0].alt   =  0;
  gcells[0].top   =  0;
  gcells[0].bot   =  0;
  gcells[0].lhs   = -1;
  gcells[0].rhs   = -1;
  gcells[0].left  =  0;
  gcells[0].above =  0;
  gcells[0].next  =  0;

  for (i = 1; i < n + 1; i++) {

      read_integer(inputFile, gcells[i].n);
      gcells[i].alt = (coor *) malloc(gcells[i].n * sizeof(coor));

      for (j = 0; j < gcells[i].n; j++) {
          read_integer(inputFile, gcells[i].alt[j][0]);
          read_integer(inputFile, gcells[i].alt[j][1]);
      }

      read_integer(inputFile, gcells[i].left);
      read_integer(inputFile, gcells[i].above);
      read_integer(inputFile, gcells[i].next);
      }

  if (!feof(inputFile)) {
      read_integer(inputFile, solution);
  }
}


static void write_outputs() {
  int i, j;

    std::cout << "Minimum area = " << MIN_AREA.load() << std::endl;

    for (i = 0; i < MIN_FOOTPRINT[0]; i++) {
      for (j = 0; j < MIN_FOOTPRINT[1]; j++) {
          if (BEST_BOARD[i][j] == 0) {bots_message(" ");}
          else                       bots_message("%c", 'A' + BEST_BOARD[i][j] - 1);
      }
      bots_message("\n");
    }
}


#ifdef BLAZE_VERSION

typedef uab::Void result_type;
// typedef size_t result_type;

struct floorplan_task
{
  enum { AREA_TASK = -1 };

  floorplan_task(int num, coor footprn, ibrd brd, int x, int nws_j0, int nws_j1, cell_ptr cll)
  : id(num), i(x), NWS_j0(nws_j0), NWS_j1(nws_j1), CELLS(cll)
  {
    memcpy(FOOTPRINT, footprn, sizeof(coor));
    memcpy(BOARD, brd, sizeof(ibrd));
  }

  floorplan_task(int num, coor footprn, ibrd brd, cell_ptr cll)
  : floorplan_task(num, footprn, brd, AREA_TASK, 0, 0, cll)
  {}

  floorplan_task() {}

#if 0 /* MANUAL_AUTO_OPS */
  floorplan_task(const floorplan_task& o)
  : id(o.id), i(o.i), NWS_j0(o.NWS_j0), NWS_j1(o.NWS_j1), CELLS(o.CELLS)
  {
    memcpy(FOOTPRINT, o.FOOTPRINT, sizeof(coor));
    memcpy(BOARD,     o.BOARD,     sizeof(ibrd));
  }

  floorplan_task(floorplan_task&& o)
  : id(o.id), i(o.i), NWS_j0(o.NWS_j0), NWS_j1(o.NWS_j1), CELLS(std::move(o.CELLS))
  {
    memcpy(FOOTPRINT, o.FOOTPRINT, sizeof(coor));
    memcpy(BOARD,     o.BOARD,     sizeof(ibrd));
  }

  floorplan_task& operator=(const floorplan_task& other)
  {
    floorplan_task  tmp(std::move(*this));

    new (this) floorplan_task (other);
    return *this;
  }

  floorplan_task& operator=(floorplan_task&& other)
  {
    floorplan_task tmp(std::move(*this));

    new (this) floorplan_task (other);
    return *this;
  }
#endif /* MANUAL_AUTO_OPS */

  int      id;

  int      i;
  int      NWS_j0;
  int      NWS_j1;

  coor     FOOTPRINT;
  ibrd     BOARD;
  cell_ptr CELLS;
};

static inline
bool is_compute_task(const floorplan_task& t)
{
  return t.i != floorplan_task::AREA_TASK;
}

static std::mutex updlock;

template <class Pool>
static
result_type add_cell_task(Pool& p, floorplan_task task)
{
  /* for each possible shape */
  for (int i = 0; i < task.CELLS.get()[task.id].n; i++)
  {
    /* compute all possible locations for nw corner */
    coor NWS[DMAX];
    int  nn = starts(task.id, i, NWS, task.CELLS.get());

    /* for all possible locations */
    for (int j = 0; j < nn; j++)
    {
      task.i = i;
      task.NWS_j0 = NWS[j][0];
      task.NWS_j1 = NWS[j][1];
      p.enq(task);
    }
  }

  return result_type(1);
}

template <class Pool>
static
auto compute_task(Pool& p, floorplan_task t) -> result_type
{
  cell*                    cells = new cell[N+1];
  cell_ptr cellptr(cells);
  ibrd                     board;
  coor                     footprint;

  memcpy(cells, t.CELLS.get(), sizeof(cell)*(N+1));

  /* extent of shape */
  cells[t.id].top = t.NWS_j0;
  cells[t.id].bot = cells[t.id].top + cells[t.id].alt[t.i][0] - 1;
  cells[t.id].lhs = t.NWS_j1;
  cells[t.id].rhs = cells[t.id].lhs + cells[t.id].alt[t.i][1] - 1;

  memcpy(board, t.BOARD, sizeof(ibrd));

  /* if the cell cannot be layed down, prune search */
  if (! lay_down(t.id, board, cells))
  {
    bots_debug("Chip %d, shape %d does not fit\n", t.id, t.i);
    return result_type(1);
  }

  /* calculate new footprint of board and area of footprint */
  footprint[0] = max(t.FOOTPRINT[0], cells[t.id].bot+1);
  footprint[1] = max(t.FOOTPRINT[1], cells[t.id].rhs+1);

  int area     = footprint[0] * footprint[1];

  /* if last cell */
  if (cells[t.id].next == 0)
  {
    /* if area is minimum, update global values */
    if (area < MIN_AREA.load(std::memory_order_relaxed))
    {
      std::lock_guard<std::mutex> guard(updlock);

      if (area < MIN_AREA.load(std::memory_order_relaxed))
      {
        MIN_FOOTPRINT[0] = footprint[0];
        MIN_FOOTPRINT[1] = footprint[1];
        memcpy(BEST_BOARD, board, sizeof(ibrd));
        bots_debug("N  %d\n", MIN_AREA.load(std::memory_order_relaxed));
        MIN_AREA.store(area, std::memory_order_relaxed);
      }
    }
  }
  else if (area < MIN_AREA.load(std::memory_order_relaxed))
  {
    /* if area is less than best area */
    add_cell_task(p, floorplan_task(cells[t.id].next, footprint, board, cellptr));
  }
  else
  {
    /* if area is greater than or equal to best area, prune search */
    bots_debug("T  %d, %d\n", area, MIN_AREA.load(std::memory_order_relaxed));
  }

  return result_type(1);
}

struct TaskHandler
{
  template <class Pool>
  result_type operator()(Pool& p, floorplan_task task)
  {
    if (is_compute_task(task))
    {
      return compute_task(p, std::move(task));
    }

    return add_cell_task(p, std::move(task));
  }
};

static inline
void print(size_t res)       { std::cerr << res << std::endl; }

static inline
void print(const uab::Void&) { }

static
void add_cell_start(int id, coor FOOTPRINT, ibrd BOARD, cell* CELLS)
{
  cell_ptr ptr(CELLS);

  result_type res = uab::execute_tasks( NUMTHREADS,
                      TaskHandler(),
                      floorplan_task( id,
                                      FOOTPRINT,
                                      BOARD,
                                      ptr
                                    )
                    );

  print(res);
}

#endif /* BLAZE_TASK */

#ifdef OMP_VERSION

struct floorplan_task
{
  enum { AREA_TASK = -1 };

  floorplan_task(int num, coor footprn, ibrd brd, int x, int nws_j0, int nws_j1, cell_ptr cll)
  : id(num), i(x), NWS_j0(nws_j0), NWS_j1(nws_j1), CELLS(cll)
  {
    memcpy(FOOTPRINT, footprn, sizeof(coor));
    memcpy(BOARD, brd, sizeof(ibrd));
  }

  floorplan_task(int num, coor footprn, ibrd brd, cell_ptr cll)
  : floorplan_task(num, footprn, brd, AREA_TASK, 0, 0, cll)
  {}

  floorplan_task() {}

  int      id;

  int      i;
  int      NWS_j0;
  int      NWS_j1;

  coor     FOOTPRINT;
  ibrd     BOARD;
  cell_ptr CELLS;
};

static std::mutex updlock;

static
void add_cell_task(floorplan_task t);

static
void compute_task(floorplan_task t)
{
  cell*                    cells = new cell[N+1];
  cell_ptr cellptr(cells);
  ibrd                     board;
  coor                     footprint;

  memcpy(cells, t.CELLS.get(), sizeof(cell)*(N+1));

  /* extent of shape */
  cells[t.id].top = t.NWS_j0;
  cells[t.id].bot = cells[t.id].top + cells[t.id].alt[t.i][0] - 1;
  cells[t.id].lhs = t.NWS_j1;
  cells[t.id].rhs = cells[t.id].lhs + cells[t.id].alt[t.i][1] - 1;

  memcpy(board, t.BOARD, sizeof(ibrd));

  /* if the cell cannot be layed down, prune search */
  if (! lay_down(t.id, board, cells))
  {
    bots_debug("Chip %d, shape %d does not fit\n", t.id, t.i);
    return;
  }

  /* calculate new footprint of board and area of footprint */
  footprint[0] = max(t.FOOTPRINT[0], cells[t.id].bot+1);
  footprint[1] = max(t.FOOTPRINT[1], cells[t.id].rhs+1);

  int area     = footprint[0] * footprint[1];

  /* if last cell */
  if (cells[t.id].next == 0)
  {
    /* if area is minimum, update global values */
    if (area < MIN_AREA.load(std::memory_order_relaxed))
    {
      std::lock_guard<std::mutex> guard(updlock);

      if (area < MIN_AREA.load(std::memory_order_relaxed))
      {
        MIN_FOOTPRINT[0] = footprint[0];
        MIN_FOOTPRINT[1] = footprint[1];
        memcpy(BEST_BOARD, board, sizeof(ibrd));
        bots_debug("N  %d\n", MIN_AREA.load(std::memory_order_relaxed));
        MIN_AREA.store(area, std::memory_order_relaxed);
      }
    }
    /* if area is less than best area */
  }
  else if (area < MIN_AREA.load(std::memory_order_relaxed))
  {
    add_cell_task(floorplan_task(cells[t.id].next, footprint, board, cellptr));
    /* if area is greater than or equal to best area, prune search */
  }
  else
  {
    bots_debug("T  %d, %d\n", area, MIN_AREA.load(std::memory_order_relaxed));
  }
}

// static
void add_cell_task(floorplan_task task)
{
  /* for each possible shape */
  for (int i = 0; i < task.CELLS.get()[task.id].n; i++)
  {
    /* compute all possible locations for nw corner */
    coor NWS[DMAX];
    int  nn = starts(task.id, i, NWS, task.CELLS.get());

    /* for all possible locations */
    for (int j = 0; j < nn; j++)
    {
      task.i = i;
      task.NWS_j0 = NWS[j][0];
      task.NWS_j1 = NWS[j][1];

      #pragma omp task firstprivate(task)
      compute_task(task);
    }
  }
}

static
void add_cell_start(int id, coor FOOTPRINT, ibrd BOARD, cell* CELLS)
{
  omp_set_num_threads(NUMTHREADS);

  cell_ptr ptr(CELLS);

  #pragma omp parallel firstprivate(id, FOOTPRINT, BOARD, ptr)
  #pragma omp single
  {
    #pragma omp taskgroup
    {
      add_cell_task(floorplan_task(id, FOOTPRINT, BOARD, ptr));
    }
  }
}

#endif /* OMP_VERSION */

#ifdef TBB_VERSION

struct floorplan_task
{
  enum { AREA_TASK = -1 };

  floorplan_task(int num, coor footprn, ibrd brd, int x, int nws_j0, int nws_j1, cell_ptr cll)
  : id(num), i(x), NWS_j0(nws_j0), NWS_j1(nws_j1), CELLS(cll)
  {
    memcpy(FOOTPRINT, footprn, sizeof(coor));
    memcpy(BOARD, brd, sizeof(ibrd));
  }

  floorplan_task(int num, coor footprn, ibrd brd, cell_ptr cll)
  : floorplan_task(num, footprn, brd, AREA_TASK, 0, 0, cll)
  {}

  floorplan_task() {}

  int      id;

  int      i;
  int      NWS_j0;
  int      NWS_j1;

  coor     FOOTPRINT;
  ibrd     BOARD;
  cell_ptr CELLS;
};

static std::mutex updlock;

static
void add_cell_task(tbb::task_group& g, floorplan_task t);

static
void compute_task(tbb::task_group& g, floorplan_task t)
{
  cell*                    cells = new cell[N+1];
  cell_ptr cellptr(cells);
  ibrd                     board;
  coor                     footprint;

  memcpy(cells, t.CELLS.get(), sizeof(cell)*(N+1));

  /* extent of shape */
  cells[t.id].top = t.NWS_j0;
  cells[t.id].bot = cells[t.id].top + cells[t.id].alt[t.i][0] - 1;
  cells[t.id].lhs = t.NWS_j1;
  cells[t.id].rhs = cells[t.id].lhs + cells[t.id].alt[t.i][1] - 1;

  memcpy(board, t.BOARD, sizeof(ibrd));

  /* if the cell cannot be layed down, prune search */
  if (! lay_down(t.id, board, cells))
  {
    bots_debug("Chip %d, shape %d does not fit\n", t.id, t.i);
    return;
  }

  /* calculate new footprint of board and area of footprint */
  footprint[0] = max(t.FOOTPRINT[0], cells[t.id].bot+1);
  footprint[1] = max(t.FOOTPRINT[1], cells[t.id].rhs+1);

  int area     = footprint[0] * footprint[1];

  /* if last cell */
  if (cells[t.id].next == 0)
  {
    /* if area is minimum, update global values */
    if (area < MIN_AREA.load(std::memory_order_relaxed))
    {
      std::lock_guard<std::mutex> guard(updlock);

      if (area < MIN_AREA.load(std::memory_order_relaxed))
      {
        MIN_FOOTPRINT[0] = footprint[0];
        MIN_FOOTPRINT[1] = footprint[1];
        memcpy(BEST_BOARD, board, sizeof(ibrd));
        bots_debug("N  %d\n", MIN_AREA.load(std::memory_order_relaxed));
        MIN_AREA.store(area, std::memory_order_relaxed);
      }
    }
    /* if area is less than best area */
  }
  else if (area < MIN_AREA.load(std::memory_order_relaxed))
  {
    add_cell_task(g, floorplan_task(cells[t.id].next, footprint, board, cellptr));
    /* if area is greater than or equal to best area, prune search */
  }
  else
  {
    bots_debug("T  %d, %d\n", area, MIN_AREA.load(std::memory_order_relaxed));
  }
}

// static
void add_cell_task(tbb::task_group& g, floorplan_task task)
{
  /* for each possible shape */
  for (int i = 0; i < task.CELLS.get()[task.id].n; i++)
  {
    /* compute all possible locations for nw corner */
    coor NWS[DMAX];
    int  nn = starts(task.id, i, NWS, task.CELLS.get());

    /* for all possible locations */
    for (int j = 0; j < nn; j++)
    {
      task.i = i;
      task.NWS_j0 = NWS[j][0];
      task.NWS_j1 = NWS[j][1];

      g.run([&g, task]() -> void { compute_task(g, task); });
    }
  }
}


static
void add_cell_start(int id, coor FOOTPRINT, ibrd BOARD, cell* CELLS)
{
  cell_ptr                 ptr(CELLS);
  tbb::task_scheduler_init init(NUMTHREADS);
  tbb::task_group          g;

  add_cell_task(g, floorplan_task(id, FOOTPRINT, BOARD, ptr));
  g.wait();
}

#endif /* TBB_VERSION */

#ifdef CILK_VERSION

struct floorplan_task
{
  enum { AREA_TASK = -1 };

  floorplan_task(int num, coor footprn, ibrd brd, cell* cll, int x, int nws_j0, int nws_j1)
  : id(num), CELLS(cll), i(x), NWS_j0(nws_j0), NWS_j1(nws_j1)
  {
    memcpy(FOOTPRINT, footprn, sizeof(coor));
    memcpy(BOARD, brd, sizeof(ibrd));
  }

  floorplan_task(int num, coor footprn, ibrd brd, cell* cll)
  : floorplan_task(num, footprn, brd, cll, AREA_TASK, 0, 0)
  {}

  floorplan_task() {}

  int   id;
  cell* CELLS;

  int   i;
  int   NWS_j0;
  int   NWS_j1;

  coor  FOOTPRINT;
  ibrd  BOARD;
};

static std::mutex updlock;

static
void add_cell_task(floorplan_task t);

// \note count tasks
// static
// void add_cell_task(floorplan_task t, cilk::reducer_opadd<size_t>& sum);
// 
// static
// void compute_task(floorplan_task t, cilk::reducer_opadd<size_t>& sum)


static
void compute_task(floorplan_task t)
{
  cell* cells = new cell[N+1];
  ibrd  board;
  coor  footprint;

  memcpy(cells, t.CELLS, sizeof(cell)*(N+1));

  /* extent of shape */
  cells[t.id].top = t.NWS_j0;
  cells[t.id].bot = cells[t.id].top + cells[t.id].alt[t.i][0] - 1;
  cells[t.id].lhs = t.NWS_j1;
  cells[t.id].rhs = cells[t.id].lhs + cells[t.id].alt[t.i][1] - 1;

  memcpy(board, t.BOARD, sizeof(ibrd));

  /* if the cell cannot be layed down, prune search */
  if (! lay_down(t.id, board, cells))
  {
    bots_debug("Chip %d, shape %d does not fit\n", t.id, t.i);
    // sum+=1;
    delete[] cells;
    return;
  }

  /* calculate new footprint of board and area of footprint */
  footprint[0] = max(t.FOOTPRINT[0], cells[t.id].bot+1);
  footprint[1] = max(t.FOOTPRINT[1], cells[t.id].rhs+1);

  int area     = footprint[0] * footprint[1];

  /* if last cell */
  if (cells[t.id].next == 0)
  {
    /* if area is minimum, update global values */
    if (area < MIN_AREA.load(std::memory_order_relaxed))
    {
      std::lock_guard<std::mutex> guard(updlock);

      if (area < MIN_AREA.load(std::memory_order_relaxed))
      {
        MIN_FOOTPRINT[0] = footprint[0];
        MIN_FOOTPRINT[1] = footprint[1];
        memcpy(BEST_BOARD, board, sizeof(ibrd));
        bots_debug("N  %d\n", MIN_AREA.load(std::memory_order_relaxed));
        MIN_AREA.store(area, std::memory_order_relaxed);
      }
    }
    /* if area is less than best area */
  }
  else if (area < MIN_AREA.load(std::memory_order_relaxed))
  {
    add_cell_task(floorplan_task(cells[t.id].next, footprint, board, cells));
    // add_cell_task(floorplan_task(cells[t.id].next, footprint, board, cells), sum);
    /* if area is greater than or equal to best area, prune search */
  }
  else
  {
    bots_debug("T  %d, %d\n", area, MIN_AREA.load(std::memory_order_relaxed));
  }

  // sum+=1;
  delete[] cells;
}

// static
// void add_cell_task(floorplan_task task, cilk::reducer_opadd<size_t>& sum)
void add_cell_task(floorplan_task task)
{
  /* for each possible shape */
  for (int i = 0; i < task.CELLS[task.id].n; i++)
  {
    /* compute all possible locations for nw corner */
    coor NWS[DMAX];
    int  nn = starts(task.id, i, NWS, task.CELLS);

    /* for all possible locations */
    for (int j = 0; j < nn; j++)
    {
      task.i = i;
      task.NWS_j0 = NWS[j][0];
      task.NWS_j1 = NWS[j][1];

      cilk_spawn compute_task(task);
      // cilk_spawn compute_task(task, sum);
    }
  }

  // sum+=1;
}

static
void set_cilk_workers(int n)
{
  assert(n <= 9999);

  char str[5];

  sprintf(str, "%d", n);

  bool success = __cilkrts_set_param("nworkers", str) == 0;
  assert(success);
}

static
void add_cell_start(int id, coor FOOTPRINT, ibrd BOARD, cell* CELLS)
{
  set_cilk_workers(NUMTHREADS);
  //   cilk::reducer_opadd<size_t> sum;

  add_cell_task(floorplan_task(id, FOOTPRINT, BOARD, CELLS));
  // add_cell_task(floorplan_task(id, FOOTPRINT, BOARD, CELLS), sum);
  // std::cerr << sum.get_value() << std::endl;
}

#endif /* CILK_VERSION */

#if QTHREADS_VERSION

template <class V>
void setenv(std::string name, V val)
{
  std::stringstream str;

  str << name << "=" << val;

  char* envset = new char[str.str().size()+1];

  memcpy(envset, str.str().c_str(), str.str().size()+1);
  putenv(envset);

  delete[] envset;
}

struct floorplan_task
{
  enum { AREA_TASK = -1 };

  floorplan_task(int num, coor footprn, ibrd brd, int x, int nws_j0, int nws_j1, qt_sinc_t* snc, cell_ptr cll)
  : id(num), i(x), NWS_j0(nws_j0), NWS_j1(nws_j1), sinc(snc), CELLS(cll)
  {
    memcpy(FOOTPRINT, footprn, sizeof(coor));
    memcpy(BOARD, brd, sizeof(ibrd));
  }

  floorplan_task(int num, coor footprn, ibrd brd, qt_sinc_t* snc, cell_ptr cll)
  : floorplan_task(num, footprn, brd, AREA_TASK, 0, 0, snc, cll)
  {}

  floorplan_task() {}

  int        id;
  int        i;
  int        NWS_j0;
  int        NWS_j1;

  coor       FOOTPRINT;
  ibrd       BOARD;

  qt_sinc_t* sinc;
  cell_ptr   CELLS;
};

static std::mutex updlock;

static
void add_cell_task(floorplan_task task);

static
aligned_t compute_task(void* qtsk)
{
  floorplan_task&          t = *reinterpret_cast<floorplan_task*>(qtsk);
  cell*                    cells = new cell[N+1];
  cell_ptr          cellptr(cells);
  ibrd                     board;
  coor                     footprint;

  memcpy(cells, t.CELLS.get(), sizeof(cell)*(N+1));

  /* extent of shape */
  cells[t.id].top = t.NWS_j0;
  cells[t.id].bot = cells[t.id].top + cells[t.id].alt[t.i][0] - 1;
  cells[t.id].lhs = t.NWS_j1;
  cells[t.id].rhs = cells[t.id].lhs + cells[t.id].alt[t.i][1] - 1;

  memcpy(board, t.BOARD, sizeof(ibrd));

  /* if the cell cannot be layed down, prune search */
  if (! lay_down(t.id, board, cells))
  {
    bots_debug("Chip %d, shape %d does not fit\n", t.id, t.i);
    qt_sinc_submit(t.sinc, nullptr);
    delete &t;
    return aligned_t();
  }

  /* calculate new footprint of board and area of footprint */
  footprint[0] = max(t.FOOTPRINT[0], cells[t.id].bot+1);
  footprint[1] = max(t.FOOTPRINT[1], cells[t.id].rhs+1);

  int area     = footprint[0] * footprint[1];

  /* if last cell */
  if (cells[t.id].next == 0)
  {
    /* if area is minimum, update global values */
    if (area < MIN_AREA.load(std::memory_order_relaxed))
    {
      std::lock_guard<std::mutex> guard(updlock);

      if (area < MIN_AREA.load(std::memory_order_relaxed))
      {
        MIN_FOOTPRINT[0] = footprint[0];
        MIN_FOOTPRINT[1] = footprint[1];
        memcpy(BEST_BOARD, board, sizeof(ibrd));
        bots_debug("N  %d\n", MIN_AREA.load(std::memory_order_relaxed));
        MIN_AREA.store(area, std::memory_order_relaxed);
      }
    }
    /* if area is less than best area */
  }
  else if (area < MIN_AREA.load(std::memory_order_relaxed))
  {
    add_cell_task(floorplan_task(cells[t.id].next, footprint, board, t.sinc, cellptr));
    /* if area is greater than or equal to best area, prune search */
  }
  else
  {
    bots_debug("T  %d, %d\n", area, MIN_AREA.load(std::memory_order_relaxed));
  }

  qt_sinc_submit(t.sinc, nullptr);
  delete &t;
  return aligned_t();
}

// static
void add_cell_task(floorplan_task task)
{
  /* for each possible shape */
  for (int i = 0; i < task.CELLS.get()[task.id].n; i++)
  {
    /* compute all possible locations for nw corner */
    coor NWS[DMAX];
    int  nn = starts(task.id, i, NWS, task.CELLS.get());

    /* for all possible locations */
    qt_sinc_expect(task.sinc, nn);
    for (int j = 0; j < nn; j++)
    {
      task.i = i;
      task.NWS_j0 = NWS[j][0];
      task.NWS_j1 = NWS[j][1];

      qthread_fork(compute_task, new floorplan_task(task), nullptr);
    }
  }
}


static
void add_cell_start(int id, coor FOOTPRINT, ibrd BOARD, cell* CELLS)
{
  cell_ptr    ptr(CELLS);
  qt_sinc_t*               sinc   = qt_sinc_create(0, nullptr, nullptr, 0);

  add_cell_task(floorplan_task(id, FOOTPRINT, BOARD, sinc, ptr));
  qt_sinc_wait(sinc, nullptr);
}


#endif /* QTHREADS_VERSION */

#ifdef BOTS_VERSION

// \note original code returns the number of explored solutions; since we do not
//       compute those in the task based model, we also do not compute the
//       compute those BOTS version.
static
void add_cell(int id, coor FOOTPRINT, ibrd BOARD, struct cell *CELLS)
{
  int  i, j, nn, area;

  ibrd board;
  coor footprint, NWS[DMAX];

  /* for each possible shape */
  for (i = 0; i < CELLS[id].n; i++) {
    /* compute all possible locations for nw corner */
    nn = starts(id, i, NWS, CELLS);
    /* for all possible locations */
    for (j = 0; j < nn; j++) {
      #pragma omp task untied private(board, footprint,area) \
      firstprivate(NWS,i,j,id,nn) \
      shared(FOOTPRINT,BOARD,CELLS,MIN_AREA,MIN_FOOTPRINT,N,BEST_BOARD)
      {
        // \was C but not C++  -> cell cells[N+1];
        cell* cells = new cell[N+1];
        memcpy(cells,CELLS,sizeof(struct cell)*(N+1));
        /* extent of shape */
        cells[id].top = NWS[j][0];
        cells[id].bot = cells[id].top + cells[id].alt[i][0] - 1;
        cells[id].lhs = NWS[j][1];
        cells[id].rhs = cells[id].lhs + cells[id].alt[i][1] - 1;

        memcpy(board, BOARD, sizeof(ibrd));

        /* if the cell cannot be layed down, prune search */
        if (! lay_down(id, board, cells))
        {
          bots_debug("Chip %d, shape %d does not fit\n", id, i);
          goto _end;
        }

        /* calculate new footprint of board and area of footprint */
        footprint[0] = max(FOOTPRINT[0], cells[id].bot+1);
        footprint[1] = max(FOOTPRINT[1], cells[id].rhs+1);
        area         = footprint[0] * footprint[1];

        /* if last cell */
        if (cells[id].next == 0)
        {
          /* if area is minimum, update global values */
          if (area < MIN_AREA) {
            #pragma omp critical
            if (area < MIN_AREA) {
              MIN_AREA         = area;
              MIN_FOOTPRINT[0] = footprint[0];
              MIN_FOOTPRINT[1] = footprint[1];
              memcpy(BEST_BOARD, board, sizeof(ibrd));
              bots_debug("N  %d\n", MIN_AREA.load(std::memory_order_relaxed));
            }
          }

          /* if area is less than best area */
        } else if (area < MIN_AREA)
        {
          add_cell(cells[id].next, footprint, board,cells);
          /* if area is greater than or equal to best area, prune search */
        } else {
          bots_debug("T  %d, %d\n", area, MIN_AREA.load(std::memory_order_relaxed));
        }
_end:;
        delete[] cells;
      }
    }
  }
  #pragma omp taskwait
}


static
void add_cell_start(int id, coor FOOTPRINT, ibrd BOARD, struct cell *CELLS)
{
    omp_set_num_threads(NUMTHREADS);

    #pragma omp parallel
    #pragma omp single
    #pragma omp taskgroup
    {
      add_cell(id, FOOTPRINT, BOARD, CELLS);
    }
}

#endif /* BOTS_VERSION */



ibrd board;

void floorplan_init (const char *filename)
{
    int i,j;

    inputFile = fopen(filename, "r");

    if(NULL == inputFile) {
        std::cerr << "Couldn't open " << filename << "file for reading" << std::endl;
        exit(1);
    }

    /* read input file and initialize global minimum area */
    read_inputs();
    MIN_AREA = ROWS * COLS;

    /* initialize board is empty */
    for (i = 0; i < ROWS; i++)
    for (j = 0; j < COLS; j++) board[i][j] = 0;
}

void compute_floorplan (void)
{
    coor footprint;
    /* footprint of initial board is zero */
    footprint[0] = 0;
    footprint[1] = 0;
    bots_message("Computing floorplan ");

    add_cell_start(1, footprint, board, gcells);

    bots_message(" completed!\n");
}

void floorplan_end (void)
{
    /* write results */
    write_outputs();
}

int floorplan_verify (void)
{
    if (solution != -1 )
      return MIN_AREA == solution ? BOTS_RESULT_SUCCESSFUL : BOTS_RESULT_UNSUCCESSFUL;
    else
      return BOTS_RESULT_NA;
}

int main(int argc, char** argv)
{
  typedef std::chrono::time_point<std::chrono::system_clock> time_point;

  std::string inputfile = "data/input.15";

  if (argc > 1) inputfile = argv[1];

  std::cout << "loading " << inputfile << std::endl;

  floorplan_init(inputfile.c_str());

#if QTHREADS_VERSION
  init_qthreads(NUMTHREADS, 20000);
#endif /* QTHREADS_VERSION */

  time_point     starttime = std::chrono::system_clock::now();
  compute_floorplan();
  time_point     endtime = std::chrono::system_clock::now();
  int            elapsedtime = std::chrono::duration_cast<std::chrono::milliseconds>(endtime-starttime).count();

  floorplan_end();

  // result validation
  int val = floorplan_verify();

  if (val != BOTS_RESULT_SUCCESSFUL)
  {
    std::cout << "result mismatch" << std::endl;
    std::cerr << "result mismatch" << std::endl;
    return 1;
  }

  // print timing
  std::cout << "time = " << elapsedtime << "ms" << std::endl;
  std::cerr << elapsedtime << std::endl;
  return 0;
}
