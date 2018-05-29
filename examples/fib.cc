// openmp:     g++ -std=c++11 -Wall -Wextra -pedantic -O2 -fopenmp -DOMP_VERSION=1 pi.cc -o /tmp/pi.bin
// blaze:      g++ -std=c++11 -Wall -Wextra -pedantic -O2 -pthread -DBLAZE_VERSION=1 pi.cc -o /tmp/pi.bin
// tbb:        g++ -std=c++11 -Wall -Wextra -pedantic -O2 -pthread -DTBB_VERSION=1 -I$TBB_HOME/include -L$TBB_HOME/lib -ltbb pi.cc -o /tmp/pi.bin
// sequential? g++ -std=c++11 -Wall -Wextra -pedantic -O2 -pthread pi.cc -o /tmp/pi.bin

#include <iostream>
#include <iomanip>
#include <chrono>
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <thread>
#include <list>
#include <vector>

#define PRINT_STATS 0

#if OMP_VERSION
#include <omp.h>
#endif

#if TBB_VERSION
#include <tbb/task_scheduler_init.h>
#include <tbb/task_group.h>
#endif

#if BLAZE_VERSION
#include "tasks.hpp"
#endif

#include "atomicutil.hpp"

#ifndef NUMTHREADS
#define NUMTHREADS (20)
#endif /* NUMTHREADS */

#ifndef FIBNUM
#define FIBNUM (6)
#endif /* NUMTHREADS */


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


template <class I>
struct fibonacci_task
{
  I num;
};


#if OMP_VERSION

//~ histogram<long double> hist(0.0, 1.0, 100);

static inline
size_t thread_num()
{
  return omp_get_thread_num();
}

template <class I>
void fib_task(fibonacci_task<I> task, uab::aligned_type<std::pair<I, size_t>, CACHELINESZ>* s)
{
  typedef fibonacci_task<I> fibonacci_task;

  size_t thrnum = thread_num();

  while (task.num > 1)
  {
    #pragma omp task
    fib_task(fibonacci_task{task.num-2}, s);

    task = fibonacci_task{task.num-1};
  }

  s[thrnum].val.first  += task.num;
  s[thrnum].val.second += 1;
}

template <class I>
auto fib_task(I num) -> std::pair<I, size_t>
{
  uab::aligned_type<std::pair<I, size_t>, CACHELINESZ> s[NUMTHREADS];
  I                                                    res  = 0;
  size_t                                               ctr  = 1;

  #pragma omp parallel firstprivate(num) shared(res,ctr,s)
  {
    size_t thrnum = thread_num();

    assert(thrnum < NUMTHREADS);
    s[thrnum] = std::make_pair(I(0), size_t(0));

    #pragma omp single
    {
      #pragma omp taskgroup
      {
        fib_task(fibonacci_task<I>{num}, s);
      }
    }

    #pragma omp atomic
    res += s[thrnum].val.first;

    #pragma omp atomic
    ctr += s[thrnum].val.second;
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

template <class D> struct globals
{
  static uab::aligned_atomic_type<D, CACHELINESZ>      result;
  static uab::aligned_atomic_type<size_t, CACHELINESZ> splits;

  static void add(D d)
  {
    D v = result.val.load(std::memory_order_relaxed);

    while (!result.val.compare_exchange_strong(v, v+d, std::memory_order_relaxed, std::memory_order_relaxed));
  }
};

template <class D>
uab::aligned_atomic_type<D, CACHELINESZ> globals<D>::result(0);


template <class D>
uab::aligned_atomic_type<size_t, CACHELINESZ> globals<D>::splits(1);


template <class G, class I>
auto compute_fib(G& taskgroup, I task) -> void // std::pair<D, size_t>
{
  while (task > 1)
  {
    taskgroup.run( [&taskgroup, task]()->void
                   { compute_fib(taskgroup, task-2);
                   }
                 );
    task =    task-1;
  }
  // return task;
}

template <class I>
auto fib_task(I num) -> std::pair<I,size_t>
{
  tbb::task_group g;

  compute_fib(g, num);

  g.wait();
  return std::pair<I,size_t>();
}

#endif /* TBB_VERSION */


#if BLAZE_VERSION

template <class I>
struct compute_fib
{
  template <class T, class A>
  I operator()(uab::pool<T,A>& tasks, T task)
  {
    typedef fibonacci_task<I> fibonacci_task;

    while (task.num > 1)
    {
      tasks.enq(fibonacci_task{task.num-2});
      task =    fibonacci_task{task.num-1};
    }

    return task.num;
  }
};

template <class I>
std::pair<I, size_t> fib_task(I num)
{
  typedef fibonacci_task<I> fibonacci_task;

  compute_fib<I> fun;

  return std::make_pair(uab::execute_tasks(NUMTHREADS, fun, fibonacci_task{ num }), 0);
}

#endif /* BLAZE_VERSION */

int main()
{
  typedef std::chrono::time_point<std::chrono::system_clock> time_point;

  time_point     starttime = std::chrono::system_clock::now();

  // executes loop in parallel
  //   and uses a reduction algorithm to combine all pi values that were
  //   computed across threads.
  #if OMP_VERSION
  omp_set_num_threads(NUMTHREADS);

  std::pair<fib_type, size_t> num  = fib_task(FIBNUM);
  #endif

  #if TBB_VERSION
  tbb::task_scheduler_init init(NUMTHREADS);

  std::pair<fib_type, size_t> num  = fib_task(FIBNUM);
  #endif

  #if BLAZE_VERSION
  std::pair<fib_type, size_t> num  = fib_task(FIBNUM);
  #endif

  time_point     endtime = std::chrono::system_clock::now();
  int            elapsedtime = std::chrono::duration_cast<std::chrono::milliseconds>(endtime-starttime).count();

  std::cout << "num    = " << FIBNUM      << std::endl;
  std::cout << "fib    = " << num.first   << std::endl;
  std::cout << "leaves = " << num.second  << std::endl;
  std::cout << "time   = " << elapsedtime << "ms" << std::endl;
  std::cerr << elapsedtime << std::endl;
  return 0;
}
