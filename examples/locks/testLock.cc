
#include <thread>
#include <iostream>
#include <sstream>
#include <list>
#include <vector>
#include <iomanip>
#include <chrono>
#include <mutex>
#include <cassert>

#ifndef UCL_RUNTIME_DATA
#define UCL_RUNTIME_DATA 1
#endif

#include "ucl/spinlock.hpp"

#include "lock-selection.hpp"

#ifndef NUM_ITER
#define NUM_ITER 1000000
#endif /* NUM_ITER */

#ifndef NUM_THREADS
#define NUM_THREADS 4
#endif /* NUM_THREADS */

typedef std::chrono::time_point<std::chrono::system_clock> time_point;

static default_lock        mtx;
static std::atomic<size_t> waiting_threads;
static std::atomic<size_t> active_threads;
static size_t              ctrlvar(0);

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


/// computes number of operations that a thread will carry out
static
size_t opsPerThread(size_t maxthread, size_t totalops, size_t thrid)
{
  assert(maxthread > thrid);

  size_t numops = totalops / maxthread;
  size_t remops = totalops % maxthread;

  if (thrid < remops) ++numops;

  return numops;
}

static
void sync_start()
{
  waiting_threads.fetch_sub(1, std::memory_order_relaxed);

  while (waiting_threads.load(std::memory_order_relaxed)) {}
}

void ptest(time_point* endtime, double* res, size_t* num_commits, const size_t cnt)
{
  size_t total = 0;
  sync_start();

  for (size_t i = cnt; i > 0; --i)
  {
    mtx.lock();
    total += ++ctrlvar;
    mtx.unlock();
  }

  size_t thrnum = active_threads.fetch_sub(1, std::memory_order_relaxed);

  if (thrnum == 1)
  {
    (*endtime) = std::chrono::system_clock::now();
  }

  (*res) = total / static_cast<double>(cnt);
  
#if UCL_RUNTIME_DATA
  *num_commits = ucl::num_commits;
#else
  *num_commits = 0;
#endif /* UCL_RUNTIME_DATA */ 
}

static
int duration(time_point starttime, time_point endtime)
{
  return std::chrono::duration_cast<std::chrono::milliseconds>(endtime-starttime).count();
}

void parallel_test(const size_t numthreads, const size_t numiter)
{
  std::cout << std::endl;

  std::list<std::thread> exp_threads;
  std::vector<double>    res_threads(numthreads, 0);
  std::vector<size_t>    run_data_threads(numthreads, 0);

  waiting_threads = numthreads;
  active_threads  = numthreads;

  time_point             endtime = std::chrono::system_clock::now();
  time_point             starttime = std::chrono::system_clock::now();

  // spawn
  for (size_t i = 0; i < numthreads; ++i)
  {
    size_t ops = opsPerThread(numthreads, numiter, i);

    exp_threads.emplace_back(ptest, &endtime, &res_threads.at(i), &run_data_threads.at(i), ops);
  }

  // join
  for (std::thread& thr : exp_threads) { thr.join(); }

  std::cout.precision(17);

  double total = 0;
  size_t total_commits = 0;

  // compute average increment to assess fairness
  for (size_t i = 0; i < numthreads; ++i)
  {
    total += res_threads.at(i);
    total_commits += run_data_threads.at(i);
    std::cout << "avg-inc " << i << ": " << std::fixed << res_threads.at(i) << std::endl;
  }

  double mean = total / numthreads;
  double variance = std::accumulate( res_threads.begin(),
                                     res_threads.end(),
                                     0.0,
                                     [mean](double sum, double v) -> double
                                     {
                                       return sum + (v - mean) * (v - mean);
                                     }
                                   );

  std::cout << std::endl
            << "std dev = " << std::sqrt(variance) << std::endl
            << "   time = " << duration(starttime, endtime) << " @" << ctrlvar << std::endl
            << "commits = " << total_commits << " (" << (total_commits * 100.0 / ctrlvar) << ")" 
            << std::endl;  
}

int main(int argc, char** args)
{
  size_t num_iter    = NUM_ITER;
  size_t num_threads = NUM_THREADS;
  size_t num_runs    = 1;

  if (argc > 1) num_iter    = asNum(*(args+1));
  if (argc > 2) num_threads = asNum(*(args+2));
  if (argc > 3) num_runs    = asNum(*(args+3));

  std::cout << "*** lock test " << num_iter << "<ops  thrds>" << num_threads << std::endl;
  std::cout << typeid(default_lock).name() << std::endl;

  while (num_runs > 0)
  {
    parallel_test(num_threads, num_iter);
    --num_runs;
  }
  return 0;
}
