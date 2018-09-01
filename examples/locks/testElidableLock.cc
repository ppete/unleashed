// g++ -ggdb -O2 -std=c++11 -Wall -Wextra -pedantic -pthread -mrtm -DHTM_ENABLED=1 -I../include testElidable.cc  -o tmp/testElidable.bin

#include <thread>
#include <iostream>
#include <sstream>
#include <list>
#include <iomanip>

#include "generalutil.hpp"
#include "locks.hpp"

#ifndef PNOITER
#define PNOITER 1000
#endif /* PNOITER */

#ifndef NUM_THREADS
#define NUM_THREADS 4
#endif /* NUM_THREADS */

#ifndef NUM_RUNS
#define NUM_RUNS 1
#endif /* NUM_RUNS */

#ifndef TEST_CONFLICTS
#define TEST_CONFLICTS 0
#endif /* TEST_CONFLICTS */


/// counts number of threads that are ready to run
static size_t waiting_threads;

/// waits until all threads have been created and are ready to run
static
void sync_start()
{
  assert(waiting_threads);
  __sync_fetch_and_sub(&waiting_threads, 1);

  while (waiting_threads) __sync_synchronize();
}

static size_t         total_time(0);
static size_t         counter(0);
static uab::ttas_lock lock;

static
void test()
{
  auto guard = uab::elide_guard(5, lock);

  ++counter;
}

typedef std::chrono::time_point<std::chrono::system_clock> time_point;
static time_point starttime;

void ptest(size_t numiter)
{
  sync_start();

  while (numiter)
  {
    test();
    --numiter;
  }
}


void parallel_test(const size_t cntthreads, const size_t cntoper)
{
  std::cout << std::endl;

  std::list<std::thread>  exp_threads;

  waiting_threads = cntthreads;

  starttime = std::chrono::system_clock::now();

  // spawn
  for (size_t i = 0; i < cntthreads; ++i)
  {
    exp_threads.emplace_back(ptest, cntoper);
  }

  // join
  for (std::thread& thr : exp_threads) thr.join();

  time_point     endtime = std::chrono::system_clock::now();
  int            elapsedtime = std::chrono::duration_cast<std::chrono::milliseconds>(endtime-starttime).count();

  __sync_synchronize();
  std::cout << "elapsed time = " << elapsedtime << "ms" << std::endl;
  std::cout << "counter = " << counter << std::endl;
  if (counter != cntthreads * cntoper) std::cout << "ERROR: counter not atomic!!" << std::endl;

  total_time += elapsedtime;

  std::cout << elapsedtime << std::endl;
  std::cout << "|| o&o." << std::endl;
}

template <class T>
static inline
size_t asNum(const T& t)
{
  size_t            res = 0;
  std::stringstream str;

  str << t;
  str >> res;

  return res;
}

int main(int argc, char** args)
{
  size_t pnoiter = PNOITER;
  size_t num_threads = NUM_THREADS;
  size_t num_runs = NUM_RUNS;

  if (argc > 1) pnoiter = asNum(*(args+1));
  if (argc > 2) num_threads = asNum(*(args+2));
  if (argc > 3) num_runs = asNum(*(args+3));

  // sequential_test(pnoiter);

  try
  {
    std::cout << "*** skiplist test " << pnoiter << "<ops  thrds>" << num_threads << std::endl;

    for (size_t i = 0; i < num_runs; ++i)
    {
      std::cout << "***** ***** ***** restart ***** " << i << std::endl;
      parallel_test(num_threads, pnoiter);
    }

    std::cerr <<"Average time: "<<total_time/num_runs<<std::endl;
    std::cout << std::endl;
  }
  catch(...)
  {
    std::cout << "Error caught..." << std::endl;
  }

  return 0;
}
