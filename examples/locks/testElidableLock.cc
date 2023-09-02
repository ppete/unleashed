// g++ -ggdb -O2 -std=c++11 -Wall -Wextra -pedantic -pthread -mrtm -DHTM_ENABLED=1 -I../include testElidable.cc  -o tmp/testElidable.bin

#include <iostream>
#include <sstream>
#include <list>
#include <iomanip>

#ifndef UCL_RUNTIME_DATA
#define UCL_RUNTIME_DATA 1
#endif

#include "ucl/unused.hpp"
#include "ucl/spinlock.hpp"
#include "ucl/thread.hpp"

#include "lock-selection.hpp"

#ifndef PNOITER
#define PNOITER 1000
#endif /* PNOITER */

#ifndef NUM_THREADS
#define NUM_THREADS 64
#endif /* NUM_THREADS */

#ifndef NUM_RUNS
#define NUM_RUNS 1
#endif /* NUM_RUNS */

#ifndef NUM_RETRIES
#define NUM_RETRIES 5
#endif /* NUM_RUNS */




/// counts number of threads that are ready to run
static std::atomic<size_t> waiting_threads;

/// waits until all threads have been created and are ready to run
static
void sync_start()
{
  assert(waiting_threads.load(std::memory_order_relaxed));

  waiting_threads.fetch_sub(1, std::memory_order_relaxed);

  while (waiting_threads.load(std::memory_order_relaxed));
}

static size_t         total_time(0);
static default_lock lock;

static
void test(size_t& counter, size_t retries)
{
  const auto& guard = ucl::elide_guard(retries, lock);

  ++counter;
}

void ptest(size_t numiter, size_t retries, size_t& counter, size_t& my_commits)
{
  sync_start();

  while (numiter)
  {
    test(counter, retries);
    --numiter;
  }
  
#if UCL_RUNTIME_DATA
  my_commits = ucl::num_commits;
#else
  my_commits = 0;
#endif /* UCL_RUNTIME_DATA */  
}


void parallel_test(const size_t cntthreads, const size_t num_retries, const size_t cntoper)
{
  typedef std::chrono::time_point<std::chrono::system_clock> time_point;

  waiting_threads = cntthreads;

  std::cout << std::endl;
  
  size_t                 counter = 0;
  std::list<ucl::thread> exp_threads;
  std::vector<size_t>    runtime_data(cntthreads, 0);
  time_point             starttime = std::chrono::system_clock::now();

  // spawn
  for (size_t i = 0; i < cntthreads; ++i)
  {
    size_t ops = cntoper / cntthreads;
    
    if (i < cntoper % cntthreads) ++ops;
    
    exp_threads.emplace_back(ptest, ops, num_retries, std::ref(counter), std::ref(runtime_data.at(i)));
  }

  // join
  for (ucl::thread& thr : exp_threads) thr.join();

  time_point     endtime = std::chrono::system_clock::now();
  int            elapsedtime = std::chrono::duration_cast<std::chrono::milliseconds>(endtime-starttime).count();

  size_t total_commits = std::accumulate(runtime_data.begin(), runtime_data.end(), 0, std::plus<size_t>());

  std::cout << "elapsed time = " << elapsedtime << "ms" << std::endl;
  
#if UCL_RUNTIME_DATA
  std::cout << "commits = " << total_commits << " (" << (total_commits * 100.0 / cntoper) << "%)" << std::endl;
#else
  std::cout << "commits = N/A (enable UCL_RUNTIME_DATA)" << std::endl;
#endif /* UCL_RUNTIME_DATA */
  
  std::cout << "counter = " << counter 
            << std::endl;
            
  if (counter != cntoper) std::cerr << "ERROR: counter not atomic!!" << std::endl;

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
  size_t num_retries = NUM_RETRIES;
  size_t num_runs = NUM_RUNS;

  if (argc > 1) pnoiter = asNum(*(args+1));
  if (argc > 2) num_threads = asNum(*(args+2));
  if (argc > 3) num_retries = asNum(*(args+3));
  if (argc > 4) num_runs = asNum(*(args+4));

  // sequential_test(pnoiter);

  try
  {
    std::cout << "*** elideguard test " << typeid(default_lock).name() << std::endl;
    std::cout << "* total-ops = " << pnoiter << std::endl;
    std::cout << "* threads   = " << num_threads << std::endl;
    std::cout << "* retries   = " << num_retries << std::endl;
    std::cout << "* runs      = " << num_runs << std::endl << std::endl;

    for (size_t i = 0; i < num_runs; ++i)
    {
      std::cout << "***** ***** ***** restart ***** " << i << std::endl;
      parallel_test(num_threads, num_retries, pnoiter);
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
