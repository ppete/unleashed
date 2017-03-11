// g++ -ggdb -O2 -std=c++11 -Wall -Wextra -pedantic -pthread -DTEST_NO_MANAGER=1 -DTEST_LOCKFREE_SKIPLIST=1 -DWITHOUT_GC=1 -I../include parSkipListUniqueElems.cc  -o tmp/parSkipListUniqueElems.bin

#include <thread>
#include <iostream>
#include <sstream>
#include <list>
#include <iomanip>

#include "generalutil.hpp"

#ifndef WITHOUT_GC
  #define GC_THREADS 1
  // #define GC_DEBUG 1
  #include <gc/gc.h>

  #include "gc-cxx11/gc_cxx11.hpp" // use GC allocator modified to work with C++11
#endif /* WITHOUT_GC */

#include "skiplist.hpp"

#ifndef PNOITER
#define PNOITER 1000
#endif /* PNOITER */

#ifndef NUM_THREADS
#define NUM_THREADS 4
#endif /* NUM_THREADS */

#ifndef NUM_RUNS
#define NUM_RUNS 1
#endif /* NUM_RUNS */

namespace lf = lockfree;
namespace fg = locking;

#if defined TEST_NO_MANAGER

template <class T>
using default_alloc = lf::just_alloc<T>;

#elif defined TEST_GC_MANAGER

template <class T>
using default_alloc = lf::gc_manager<T,  gc_allocator_cxx11>;

#elif defined TEST_EPOCH_MANAGER

template <class T>
using default_alloc = lf::epoch_manager<T>;

#elif defined TEST_PUB_SCAN_MANAGER

template <class T>
using default_alloc = lf::pub_scan_manager<T>;

#else /* undefined MANAGER */

  #error "preprocessor define for memory manager is needed (TEST_JUST_ALLOC, TEST_GC_MANAGER, TEST_EPOCH_MANAGER, TEST_PUB_SCAN_MANAGER)"

#endif /* TEST_NO_MANAGER */

#ifndef TEST_MAX_LEVELS
#define TEST_MAX_LEVELS 32
#endif

#ifndef TEST_CONFLICTS
#define TEST_CONFLICTS 0
#endif

static const size_t LEVELS = TEST_MAX_LEVELS;

#if defined TEST_LOCKING_SKIPLIST

template <class T>
using skiplist = fg::skiplist<T, std::less<T>, default_alloc<T>, LEVELS>;

#elif defined TEST_LOCKFREE_SKIPLIST

template <class T>
using skiplist = lf::skiplist<T, std::less<T>, default_alloc<T>, LEVELS>;

#else

 #error "preprocessor define for container class is needed (TEST_LOCKING_SKIPLIST, TEST_LOCKFREE_SKIPLIST)"

#endif


typedef skiplist<int> container_type;

int total_time = 0;

struct alignas(CACHELINESZ) ThreadInfo
{
  container_type* container;
  size_t          num;
  size_t          fail;
  size_t          succ;
  size_t          pnoiter;
  size_t          num_threads;
  size_t          num_held_back;
  // unsigned        cpunode_start;
  // unsigned        cpunode_end;
  // std::set        elemsin;
  // std::set        elemsout;

  ThreadInfo(container_type* r, size_t n, size_t cntiter, size_t cntthreads)
  : container(r), num(n), fail(0), succ(0), pnoiter(cntiter), num_threads(cntthreads),
    num_held_back(0)
    // , cpunode_start(0), cpunode_end(0)
  {
    assert(num < num_threads);
  }

  ThreadInfo()
  : ThreadInfo(nullptr, 0, 0, 0)
  {}
};

static
void fail()
{
  throw std::logic_error("error");
}

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

/// computes number of operations that a thread will carry out
static
size_t opsPerThread(size_t maxthread, size_t totalops, size_t thrid)
{
  assert(maxthread > thrid);

  size_t numops = totalops / maxthread;
  size_t remops = totalops % maxthread;

  if (remops > 0 && thrid < remops) ++numops;

  return numops;
}

/// computes number of operations in the main loop; the remainder of
///   operations will be executed prior to the main loop (i.e., insert elements)
static
size_t opsMainLoop(size_t numops)
{
  return numops - (numops / 10);
}

/// computes the expected size after all operations have been finished
///   requires deterministic skiplist operations
static
int _expectedSize(size_t maxthread, size_t totalops)
{
  int elems = 0;

  for (size_t i = 0; i < maxthread; ++i)
  {
    int numops  = opsPerThread(maxthread, totalops, i);
    int nummain = opsMainLoop(numops);

    // special case when nummain is even and small
    //   then the insert happens after the delete
    if ((numops - nummain == 0) && (nummain % 2 == 0))
      elems += numops / 2;
    else
      elems += (numops - nummain);

    // if odd, we insert one more
    if (nummain % 2) ++elems;
  }

  return elems;
}

#if TEST_CONFLICTS

//
// contention at end of list

static
size_t genElem(size_t num, size_t thrid, size_t maxthread, size_t)
{
  return num * maxthread + thrid;
}

static
int expectedSize(size_t maxthread, size_t totalops)
{
  return _expectedSize(maxthread, totalops);
}

#else

// generate elems so threads operate by and large in disjoint segments of the skiplist

static
int expectedSize(size_t maxthread, size_t totalops)
{
  return _expectedSize(maxthread, totalops);
}

static
size_t genElem(size_t num, size_t thrid, size_t, size_t)
{
  // note, the multiplier must be adjusted with higher number of elements
  return thrid * 2000000 + num;
}

#endif

std::atomic<bool> ptest_failed(false);

static
void container_test_prefix(ThreadInfo& ti)
{
  const size_t      tinum   = ti.num;
  size_t            numops  = opsPerThread(ti.num_threads, ti.pnoiter, tinum);
  const size_t      threadops = numops;
  const size_t      nummain = opsMainLoop(numops);
  size_t            wrid    = 0;

  try
  {
    while (numops > nummain)
    {
      int     elem = genElem(++wrid, tinum, ti.num_threads, threadops);
      int     succ = ti.container->insert(elem);

      assert(succ >= 0), unused(succ);
      ++ti.succ;
      // std::cout << "added " << elem << " " << succ << std::endl;

      --numops;
    }
  }
  catch (int errc)
  {
    ptest_failed = true;
    std::cout << "err: " << errc << std::endl;
  }
}


static
void container_test(ThreadInfo& ti)
{
  const size_t      tinum   = ti.num;
  size_t            numops  = opsPerThread(ti.num_threads, ti.pnoiter, tinum);
  const size_t      threadops = numops;
  const size_t      nummain = opsMainLoop(numops);
  size_t            wrid    = numops - nummain;
  size_t            rdid    = wrid / 2;

  try
  {
    // set numops to nummain (after prefix has been executed)
    numops = nummain;

    assert(numops > 0);
    while (numops)
    {
      if (numops % 2)
      {
        int    elem = genElem(++wrid, tinum, ti.num_threads, threadops);
        int    succ = ti.container->insert(elem);

        assert(succ >= 0), unused(succ);
        ++ti.succ;
      }
      else
      {
        int    elem = genElem(++rdid, tinum, ti.num_threads, threadops);
        int    succ = ti.container->erase(elem);
        if (succ > 0) ++ti.succ; else ++ti.fail;
      }

      --numops;
    }

    ti.container->getAllocator().release_memory();
    ti.num_held_back = ti.container->getAllocator().has_unreleased_memory();
  }
  catch (int errc)
  {
    ptest_failed = true;
    std::cout << "err: " << errc << std::endl;
  }
}

typedef std::chrono::time_point<std::chrono::system_clock> time_point;
static time_point starttime;

void ptest(ThreadInfo *const ti)
{
  gc_cxx_thread_context gc_guard;

  container_test_prefix(*ti);
  sync_start();

  if (ti->num == 0) starttime = std::chrono::system_clock::now();
  container_test(*ti);
}

void sequential_test(size_t cntoper)
{
  container_type cont;
  ThreadInfo     info(&cont, 0, cntoper, 1);

  std::cout << std::endl;
  container_test(info);
  std::cout << "seq OK" << std::endl;
}


void parallel_test(const size_t cntthreads, const size_t cntoper)
{
  std::cout << std::endl;

  std::list<std::thread>  exp_threads;
  std::list<ThreadInfo>   thread_info;
  container_type          cont;

  waiting_threads = cntthreads;

  // spawn
  for (size_t i = 0; i < cntthreads; ++i)
  {
    thread_info.push_back(ThreadInfo(&cont, i, cntoper, cntthreads));
    exp_threads.emplace_back(ptest, &thread_info.back());
  }

  // join
  for (std::thread& thr : exp_threads) thr.join();

  time_point     endtime = std::chrono::system_clock::now();
  int            elapsedtime = std::chrono::duration_cast<std::chrono::milliseconds>(endtime-starttime).count();

  __sync_synchronize();
  const int szlst = std::distance(cont.qbegin(), cont.qend());
  std::cout << "elapsed time = " << elapsedtime << "ms" << std::endl;
  std::cout << "skiplist size = " << szlst << std::endl;
  total_time += elapsedtime;

  std::cout << elapsedtime << std::endl;

  size_t total_succ = 0;
  size_t total_fail = 0;

  for (ThreadInfo& info : thread_info)
  {
    std::cout << "i: " << info.num << "  "
              << info.succ << "(" << (info.fail + info.succ) << ")"
              //~ << " [ " << thread_info[i].cpunode_start
              //~ << " - " << thread_info[i].cpunode_end
              << " ]  x "
              << "  ( " << info.num_held_back << " obj held )"
              << std::endl;

    total_succ += info.succ;
    total_fail += info.fail;
  }

  if (expectedSize(cntthreads, cntoper) != szlst)
  {
    std::cout << "Unexpected size " << expectedSize(cntthreads, cntoper)
              << " <exp != act> " << szlst << std::endl;
    fail();
  }

  std::cout << szlst << std::endl;
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
#ifndef WITHOUT_GC
  GC_INIT();
  GC_allow_register_threads();
  // GC_find_leak = 1;
#endif /* WITHOUT_GC */

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
    std::cout << typeid(container_type).name() << "  " << container_type::MAXLEVEL << std::endl;

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
