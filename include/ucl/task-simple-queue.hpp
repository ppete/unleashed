
#ifndef _UCL_TASK_SINGLE_QUEUE
#define _UCL_TASK_SINGLE_QUEUE 1


#include "ucl/task.hpp"
#include "ucl/queue.hpp"
#include "ucl/generalutil.hpp"

#define ENABLE_UNSAFE_LOCK_IMPLS 1
#include "ucl/locks.hpp"

namespace ucl
{
  /// simple, straight-forward single queue design
  template <class MPMCBag>
  struct single_q_pool
  {
    // define the underlying queue
    typedef MPMCBag                        bag_t;

    // defines the task type this pool is able to handle
    typedef typename bag_t::value_type     value_type;

    // allocator type
    typedef typename bag_t::allocator_type allocator_type;


    // data members
    const uint_fast32_t                    workers;
    bag_t                                  queue;
    ucl::aligned_atomic_type<int_fast32_t> activethreads;

    static thread_local bool active;

    single_q_pool(uint_fast32_t numthreads, value_type work, const allocator_type& alloc = allocator_type())
    : workers(numthreads), queue(alloc), activethreads(numthreads)
    {
      queue.enq(std::move(work));
    }

    // returns the number of worker threads
    size_t num_threads() { return workers; }

    // initializes pool thread local data for thisthread
    // and sets its parent thread
    void init(size_t parent, size_t thisthread)
    {
      unused(parent, thisthread);

      active = true;
    }

    // attempts to get a task from some queue and stores it in t
    // \result true, iff a task could be retrieved
    bool deq(value_type& t)
    {
      bool res;
      std::tie(t, res) = queue.deq();

      return res;
    }

    // adds a new task t to this thread's queue
    void enq(value_type t)
    {
      queue.enq(std::move(t));
    }

    // called before work on task starts
    void work_started()
    {
      if (!active)
      {
        activethreads.val.fetch_add(1, std::memory_order_relaxed);
        active = true;
      }
    }

    // called after work on task finished
    void work_completed() {}

    bool became_inactive()
    {
      if (!active) return false;

      active = false;
      return true;
    }

    // returs true if there is at least one active task in the system
    bool work_available()
    {
      size_t cntactv = became_inactive()
                          ? activethreads.val.fetch_sub(1, std::memory_order_relaxed) - 1
                          : activethreads.val.load(std::memory_order_relaxed)
                          ;

      return cntactv > 0;
    }

    // releases (thread local) resources after all tasks have been handled
    // \details
    //    it can be assumed that this method is called during quiescent time
    void qrelease_memory()
    {
      queue.qrelease_memory();
    }
  };

  template <class MPMCBag>
  thread_local
  bool single_q_pool<MPMCBag>::active = false;
}

#endif /* _UCL_TASK_SINGLE_QUEUE */
