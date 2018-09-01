#include <thread>
#include <iostream>
#include <sstream>
#include <list>

#include "atomicutil.hpp"
#include "list.hpp"

#ifndef PNOITER
#define PNOITER 1000
#endif /* PNOITER */

#ifndef NUM_THREADS
#define NUM_THREADS 6
#endif /* NUM_THREADS */

#ifndef NUM_RUNS
#define NUM_RUNS 1
#endif /* NUM_RUNS */

namespace lf = lockfree;

template <class T>
//~ using list = lockfree::sorted_list<T, std::less<T>, lockfree::epoc_manager<T> >;
// using list = lockfree::sorted_list<T, std::less<T>, lockfree::pub_scan_manager<T> >;
using list = lockfree::sorted_list<T, std::less<T>, lockfree::just_alloc<T> >;

typedef list<int> container_type;

struct history_entry
{
  void* ptr;
  char  opc;

  history_entry(void* p, char o)
  : ptr(p), opc(o)
  {}

  history_entry()
  : history_entry(nullptr, ' ')
  {}
};

struct alignas(CACHELINESZ) ThreadInfo
{
  static const size_t SZHISTORY = 32;

  container_type* container;
  size_t          num;
  size_t          fail;
  size_t          succ;
  size_t          pnoiter;
  size_t          num_threads;
  size_t          num_held_back;
  size_t          hcnt;
  history_entry   history[SZHISTORY];

  // unsigned        cpunode_start;
  // unsigned        cpunode_end;
  // std::set        elemsin;
  // std::set        elemsout;

  ThreadInfo(container_type* r, size_t n, size_t cntiter, size_t cntthreads)
  : container(r), num(n), fail(0), succ(0), pnoiter(cntiter), num_threads(cntthreads), hcnt(0)
  {
    assert(num < num_threads);

    // for (history_entry& el: history) el.opc = ' ';
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

static size_t waiting_threads;

static
void sync_start()
{
  assert(waiting_threads);
  __sync_fetch_and_sub(&waiting_threads, 1);

  while (waiting_threads) __sync_synchronize();
}


/// generates numbers that are more likely to cause contention
static
size_t genElem(size_t num, size_t thrid, size_t maxthread)
{
  return num * maxthread + thrid;
}

/*
/// generates numbers that are less likely to cause contention
static
size_t genElem(size_t num, size_t thrid, size_t)
{
  return thrid * 100000 + num;
}
*/

void* ptest(void* info)
{
  ThreadInfo *const ti = static_cast<ThreadInfo*>(info);
  const size_t      tinum = ti->num;
  size_t            wrid = 0;
  size_t            rdid = 0;
  size_t            totalops = ti->pnoiter / ti->num_threads;
  size_t            remops   = ti->pnoiter % ti->num_threads;

  if (remops > 0 && tinum < remops) ++totalops;

  // query node assignments
  sync_start();

  const size_t numadd = totalops - (totalops / 10);
  while (totalops > numadd)
  {
    void* succ = ti->container->insert(genElem(++wrid, tinum, ti->num_threads));
    assert(succ);

    ti->history[++ti->hcnt % ThreadInfo::SZHISTORY] = history_entry(succ, 'i');

    ++ti->succ;
    --totalops;
  }

  rdid = wrid / 2;

  assert(totalops > 0);
  while (totalops)
  {
    if (totalops % 2)
    {
      void* succ = ti->container->insert(genElem(++wrid, tinum, ti->num_threads));
      assert(succ);
      ti->history[++ti->hcnt % ThreadInfo::SZHISTORY] = history_entry(succ, 'i');
      ++ti->succ;
    }
    else
    {
      size_t elem = genElem(++rdid, tinum, ti->num_threads);
      void*  succ = ti->container->erase(elem);

      ti->history[++ti->hcnt % ThreadInfo::SZHISTORY] = history_entry(succ, 'd');
      if (succ) ++ti->succ; else ++ti->fail;
    }

    --totalops;
  }

  ti->container->get_allocator().release_memory();
  ti->num_held_back = ti->container->get_allocator().has_unreleased_memory();

  // query at end of thread
  return nullptr;
}

template <class _Cont>
void assert_order(_Cont& cont)
{
  typename _Cont::iterator   aa = cont.qbegin();
  typename _Cont::iterator   zz = cont.qend();
  typename _Cont::value_type prev = *aa;

  ++aa;
  while (aa != zz)
  {
    typename _Cont::value_type curr = *aa;

    if (prev >= curr)
    {
      std::cerr << prev << " not less than " << curr << std::endl;
    }

    prev = curr;
    ++aa;
  }
}



void parallel_test(size_t cntthreads, size_t cntoper)
{
  typedef std::chrono::time_point<std::chrono::system_clock> time_point;

  std::cout << std::endl;

  std::list<pthread_t>  exp_threads;
  std::list<ThreadInfo> thread_info;
  container_type        cont;

  waiting_threads = cntthreads;
  time_point            starttime = std::chrono::system_clock::now();

  // spawn
  for (size_t i = 0; i < cntthreads; ++i)
  {
    thread_info.push_back(ThreadInfo(&cont, i, cntoper, cntthreads));
    exp_threads.push_back(pthread_t());

    const int t_status = pthread_create(&exp_threads.back(), NULL, ptest, &thread_info.back());

    if (t_status) fail();
  }

  // join
  for (pthread_t thr : exp_threads) pthread_join(thr, NULL);

  time_point     endtime = std::chrono::system_clock::now();
  int            elapsedtime = std::chrono::duration_cast<std::chrono::milliseconds>(endtime-starttime).count();

  __sync_synchronize();
  std::cout << "elapsed time = " << elapsedtime << "ms" << std::endl;
  std::cout << "container size = " << cont.qcount() << std::endl;

  std::cerr << elapsedtime << std::endl;

  size_t total_succ = 0;
  size_t total_fail = 0;

  for (ThreadInfo& info : thread_info)
  {
    std::cout << "i: " << info.num << "  "
              << info.succ << "(" << (info.fail + info.succ) << ")"
              //~ << " [ " << thread_info[i].cpunode_start
              //~ << " - " << thread_info[i].cpunode_end
              << " ]   "
              << info.num_held_back << "mem"
              << std::endl;

    total_succ += info.succ;
    total_fail += info.fail;
  }

/*
  if (cont.get_allocator().error_state())
  {
    bool active = true;

    for (size_t sz = 1; active && sz <= ThreadInfo::SZHISTORY; ++sz)
    {
      active = false;

      for (ThreadInfo& info : thread_info)
      {
        char  opc = ' ';
        void* xxx = reinterpret_cast<void*>((int64_t)-1);

        if (info.hcnt >= sz)
        {
          opc = info.history[sz].opc;
          xxx = info.history[sz].ptr;
          active = true;
        }

        std::cerr << opc;
        std::cerr.setf(std::ios::left, std::ios::adjustfield);
        std::cerr.width(16);
        if (opc == ' ') std::cerr << ""; else std::cerr << xxx;
      }

      std::cerr << std::endl;
    }

    assert(false);
  }
  */

  assert_order(cont);

  std::cout << cont.qcount() << std::endl;
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
  long seed = time(0);

  seed = 1399724569;

  // srand(seed); // sets up the "rand()"
  std::cout << "seed = " << seed << std::endl;

  size_t pnoiter = PNOITER;
  size_t num_threads = NUM_THREADS;
  size_t num_runs = NUM_RUNS;

  if (argc > 1) pnoiter = asNum(*(args+1));
  if (argc > 2) num_threads = asNum(*(args+2));
  if (argc > 3) num_runs = asNum(*(args+3));

  try
  {
    std::cout << "*** skiplist test " << num_threads << "<thrds  ops>" << pnoiter << std::endl;
    std::cout << "    sizeof(size_t) =  " << sizeof(size_t) << std::endl;
    std::cout << "    sizeof(void*)  =  " << sizeof(void*) << std::endl;

    for (size_t i = 0; i < num_runs; ++i)
    {
      parallel_test(num_threads, pnoiter);
    }

    std::cout << std::endl;
  }
  catch(...)
  {
    std::cerr << "Error caught..." << std::endl;
  }

  return 0;
}


