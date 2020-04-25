/// \file spinlock.hpp
/// \brief File provides (various) spinlock implementation(s) in C++11
///        and methods to elide them in the presence of hardware transactional
///        memory (requires HTM_ENABLED be set to 1).
/// \todo  separate some classes that need to hold data into header and
///        implementation file.
/// \author Amalee Wilson, Peter Pirkelbauer
///
/// To seamlessly interact with systems that support transactional memory,
/// we introduce the Elidable concept. Elidable is similar to
/// std::mutex, adds is_free(), but makes try_lock optional.
/// is_free is used by elided_lock and elidable_guard to test whether the
///   execution context is transactional.
/// try_lock is optional b/c some locking techniques are unable to support it.
///
/// \code{.cpp}
/// concept Elidable
/// {
///   /// no argument constructor
///   Elidable();
///
///   /// OPTIONAL: attempts to acquire lock once
///   /// \return true if successful, false otherwise
///   bool try_lock();
///
///   /// blocks until the lock is acquired (acquire semantics)
///   void lock();
///
///   /// unlocks the lock (release semantics)
///   void unlock();
///
///   /// tests whether the lock is taken. Used to test whether code after
///   ///   lock acquisition runs in transactional context.
///   bool is_free();
/// }
/// \endcode

#ifndef _UNLEASHED_LOCKS_HPP
#define _UNLEASHED_LOCKS_HPP

#include <atomic>
#include <random>
#include <thread>
#include <limits>

#include "ucl/atomicutil.hpp"

#if HTM_ENABLED
#include "ucl/htm.hpp"
#endif /* HTM_ENABLED */

#if 0 // defined(__x86_64__)
#include "xmmintrin.h"

/// \private
static inline
void pausenow() { _mm_pause(); }

#else

/// \private
static inline
void pausenow() {}

#endif



namespace ucl
{

#if UCL_RUNTIME_DATA
  thread_local size_t num_commits  = 0;
#endif /* UCL_RUNTIME_DATA */


  /// \brief Overly simple test and set lock.
  /// \details demonstrates the problem of just using compare_and_exchange
  ///          when spinning. For a much better implementation see ttas_lock.
  struct tas_lock
  {
      typedef int_fast8_t             token_type;
      typedef std::atomic<token_type> control_type;
      typedef control_type&           native_handle_type;

      tas_lock()
      : mem(0)
      {}

      bool try_lock()
      {
        token_type oldval = 0;

        return mem.compare_exchange_weak(oldval, 1, std::memory_order_acquire, std::memory_order_relaxed);
      }

      void lock()
      {
        while (!try_lock()) ;
      }

      void unlock()
      {
        mem.store(0, std::memory_order_release);
      }

      bool is_free() const
      {
        return mem.load(std::memory_order_relaxed) == 0;
      }

      native_handle_type& native_handle()
      {
        return mem;
      }

    private:
      control_type mem;

      // disabled functions
      tas_lock(const tas_lock&) = delete;
      tas_lock(tas_lock&&) = delete;
      tas_lock& operator=(const tas_lock&) = delete;
      tas_lock& operator=(const tas_lock&&) = delete;
  };

  /// \brief ttas_lock (aka shadow lock), uses reads while spinning
  struct ttas_lock
  {
      typedef int_fast8_t             token_type;
      typedef std::atomic<token_type> control_type;
      typedef control_type&           native_handle_type;

      ttas_lock()
      : mem(0)
      {}

      bool try_lock()
      {
        token_type oldval = 0;

        return mem.compare_exchange_weak(oldval, 1, std::memory_order_acquire, std::memory_order_relaxed);
      }

      void lock()
      {
        for (;;)
        {
          while (mem.load(std::memory_order_relaxed)) { pausenow(); }

          if (try_lock()) return;
        }
      }

      void unlock()
      {
        mem.store(0, std::memory_order_release);
      }

      bool is_free() const
      {
        return (mem.load(std::memory_order_relaxed) == 0);
      }

      native_handle_type& native_handle()
      {
        return mem;
      }

    private:
      control_type mem;

      // disabled functions
      ttas_lock(const ttas_lock&) = delete;
      ttas_lock(ttas_lock&&) = delete;
      ttas_lock& operator=(const ttas_lock&) = delete;
      ttas_lock& operator=(const ttas_lock&&) = delete;
  };


#if HTM_ENABLED
  /// \brief wrapper class around a spinlock that attempts to elide the lock
  ///        and execute the locked section as transaction.
  /// \details the class does not implement try_lock (what's its meaning in a TX?).
  /// \tparam _L the type of the underlying lock
  /// \tparam NUMTRIES number of attempts a TX is started before acquiring
  ///         the lock.
  template <class _L, size_t NUMTRIES = 10>
  struct elided_lock
  {
    _L the_lock;

    elided_lock()
    : the_lock()
    {}

    /// \brief starts a transaction or acquires the lock
    /// \return true, if in TX mode; false otherwise
    // \todo better tracking of aborts so transaction length can be updated.
    bool lock()
    {
      size_t i = NUMTRIES;

      while (i > 0)
      {
        if (htm::tx::begin())
        {
          if (the_lock.is_free()) return true;

          htm::tx::abort<81>();
        }

        --i;
      }

      the_lock.lock();
      return false;
    }

    /// \brief returns true, if the lock is free
    ///        i.e., if the TX was successfully started
    bool is_free() const
    {
      return the_lock.is_free();
    }

    /// commits transaction, or unlocks the lock
    void unlock()
    {
      // if lock was not taken, TX must be active
      if (the_lock.is_free())
      {
        htm::tx::end();
#if UCL_RUNTIME_DATA
        ++num_commits;
#endif /* UCL_RUNTIME_DATA */
        return;
      }

      the_lock.unlock();
    }
  };
#endif

  /// \brief ttas-lock that backs off under contention using increasing intervals
  /// \details does not implement try_lock
  template <size_t MAX_SLEEP_TIME = (1 << 8)>
  struct ttas_lock_backoff
  {
    typedef int_fast8_t token_type;

    enum { MIN_SLEEP = 10, MAX_SLEEP = MAX_SLEEP_TIME };

    ttas_lock_backoff()
    : mem(0)
    {}

    void lock()
    {
      long       sleep_time = MIN_SLEEP;
      token_type oldval = 0;

      while (!mem.compare_exchange_weak(oldval, 1, std::memory_order_acquire, std::memory_order_relaxed))
      {
        sleep_time = sleep_time << 1;
        if (sleep_time > MAX_SLEEP) sleep_time = MAX_SLEEP;

        std::this_thread::sleep_for(std::chrono::nanoseconds(sleep_time));

        while (mem.load(std::memory_order_relaxed)) {}
        oldval = 0;
      }
    }

    bool try_lock()
    {
      if (mem.load(std::memory_order_relaxed)) return false;

      token_type oldval = 0;

      return mem.compare_exchange_weak(oldval, 1, std::memory_order_acquire, std::memory_order_relaxed);
    }

    void unlock()
    {
      mem.store(0, std::memory_order_release);
    }

    bool is_free() const
    {
      return mem.load(std::memory_order_relaxed) == 0;
    }

    private:
      std::atomic<token_type>              mem;

      // static thread_local std::minstd_rand ttas_lock_backoff_rand;

      // disabled functions
      ttas_lock_backoff(const ttas_lock_backoff&) = delete;
      ttas_lock_backoff(ttas_lock_backoff&&) = delete;
      ttas_lock_backoff& operator=(const ttas_lock_backoff&) = delete;
      ttas_lock_backoff& operator=(ttas_lock_backoff&&) = delete;
  };

  // template <size_t MAX_SLEEP_TIME>
  // thread_local std::minstd_rand
  // ttas_lock_backoff<MAX_SLEEP_TIME>::ttas_lock_backoff_rand(77);

  typedef ttas_lock_backoff<> ttas_lock_backoff_default;


  /// \private
  /// \brief    anderson ticket lock
  template <size_t ALIGNSZ=CACHELINESZ, size_t MAX_TICKETS=128>
  struct anderson_lock
  {
    typedef uint_fast16_t ticket_type;
    typedef uint64_t      counter_type;

    enum { MAX_TICKET_COUNT = MAX_TICKETS };

    anderson_lock()
    : waitctr(0), anderson_ticket(0)
    {
      static_assert( std::numeric_limits<ticket_type>::max() > MAX_TICKET_COUNT,
                     "size of ticket array is too large"
                   );

      for (size_t i = 0; i < MAX_TICKET_COUNT; ++i)
        waitlist[i].val.store(0, std::memory_order_relaxed);
    }

    void lock()
    {
      const counter_type no = waitctr.val.fetch_add(1, std::memory_order_relaxed);
      waitlist_entry&    myticket = waitlist[no % MAX_TICKET_COUNT].val;

      while (myticket.load(std::memory_order_acquire) != no) {}

      anderson_ticket = no;
    }

    void unlock()
    {
      const counter_type nextnum = (anderson_ticket+1);
      waitlist_entry&    other = waitlist[nextnum%MAX_TICKET_COUNT].val;

      other.store(nextnum, std::memory_order_release);
    }

    bool is_free() const
    {
      const counter_type no = waitctr.val.load(std::memory_order_relaxed);
      const waitlist_entry& myticket = waitlist[no % MAX_TICKET_COUNT].val;

      return myticket.load(std::memory_order_relaxed) == no;
    }

    private:
      typedef std::atomic<counter_type>             waitlist_entry;
      typedef aligned_type<waitlist_entry, ALIGNSZ> aligned_entry;

      aligned_atomic_type<counter_type>    waitctr;
      aligned_entry                        waitlist[MAX_TICKET_COUNT];
      ticket_type                          anderson_ticket;

      // disabled functions
      anderson_lock(const anderson_lock&) = delete;
      anderson_lock(anderson_lock&&) = delete;
      anderson_lock& operator=(const anderson_lock&) = delete;
      anderson_lock& operator=(anderson_lock&&) = delete;
  };


  /// \private
  struct clh_lock_node
  {
    explicit
    clh_lock_node(int val)
    : flag(val)
    {}

    std::atomic<int> flag;
  };

  /// \private
  thread_local clh_lock_node* clh_curr = nullptr;

  /// \private
  thread_local clh_lock_node* clh_next = nullptr;

  /// \private
  struct clh_lock
  {
    clh_lock()
    : next(new clh_lock_node(1)), owner_curr(nullptr)
    {}

    void lock()
    {
      // setup next object
      if (clh_next == nullptr)
        clh_next = new clh_lock_node(0);
      else
        clh_next->flag.store(0, std::memory_order_relaxed);

      // swap my node with current lock node
      clh_lock_node* clh_curr = next.val.exchange(clh_next, std::memory_order_acq_rel);

      // load until green light
      while (!clh_curr->flag.load(std::memory_order_relaxed)) {}

      // acquire fence to sync with lock acquisition
      std::atomic_thread_fence(std::memory_order_acquire);

      owner_curr = clh_curr;
    }

    void unlock()
    {
      // set free next thread
      owner_curr->flag.store(1, std::memory_order_release);

      // keep one node for the next lock acquisition
      // \todo could be improved for codes relying on nested locking
      if (clh_next) delete clh_next;

      clh_next = owner_curr;
    }

    bool is_free() const
    {
      return (  (clh_next == nullptr)
             || (clh_next->flag.load(std::memory_order_relaxed) == 1)
             );
    }

    private:
      ucl::aligned_atomic_type<clh_lock_node*> next;
      clh_lock_node*                           owner_curr;

      // disabled functions
      clh_lock(const clh_lock&) = delete;
      clh_lock(clh_lock&&) = delete;
      clh_lock& operator=(const clh_lock&) = delete;
      clh_lock& operator=(clh_lock&&) = delete;
  };

  /// \private
  struct mcs_lock_node
  {
    typedef int_fast8_t token_type;

    mcs_lock_node()
    : flag(0), next(nullptr)
    {}

    std::atomic<token_type>     flag;
    std::atomic<mcs_lock_node*> next;
  };

  /// \private
  thread_local mcs_lock_node* mcs_freelist = nullptr;

  /// \private
  struct mcs_lock
  {
    mcs_lock()
    : next(nullptr), owner(nullptr)
    {}

    mcs_lock_node* alloc_node()
    {
      if (mcs_freelist == nullptr) return new mcs_lock_node;

      mcs_lock_node* mcs_node = mcs_freelist;

      mcs_node->next.store(nullptr, std::memory_order_relaxed);
      mcs_node->flag.store(0, std::memory_order_relaxed);
      return mcs_node;
    }

    void dealloc_node(mcs_lock_node* mcs_node)
    {
      if (mcs_freelist != nullptr) delete mcs_freelist;

      mcs_freelist = mcs_node;
    }

    void lock()
    {
      mcs_lock_node* mcs_node = alloc_node();

      // \mot preceding stores must complete before we publish the node
      mcs_lock_node* prev = next.val.exchange(mcs_node, std::memory_order_release);

      if (prev != nullptr)
      {
        // \todo needed? std::atomic_thread_fence(std::memory_order_consume);
        prev->next.store(mcs_node, std::memory_order_relaxed);

        while (!mcs_node->flag.load(std::memory_order_relaxed)) {}
      }

      std::atomic_thread_fence(std::memory_order_acquire);
      owner.val = mcs_node;
    }

    void unlock()
    {
      // if our flag has been set already, we do not have to attempt c&s
      mcs_lock_node* mcs_node = owner.val;
      mcs_lock_node* old = mcs_node->next.load(std::memory_order_relaxed);

      if (old == nullptr)
      {
        old = mcs_node;
        next.val.compare_exchange_strong(old, nullptr, std::memory_order_release, std::memory_order_relaxed);
      }

      // if someone else is waiting...
      if (old != mcs_node)
      {
        // wait until the next in line has revealed itself
        while (mcs_node->next.load(std::memory_order_relaxed) == nullptr) {}

        mcs_lock_node* nxt = mcs_node->next.load(std::memory_order_consume);

        nxt->flag.store(1, std::memory_order_release);
      }
    }

    bool is_free() const
    {
      return next.val.load(std::memory_order_relaxed) == nullptr;
    }

    private:
      aligned_atomic_type<mcs_lock_node*> next;
      aligned_type<mcs_lock_node*>        owner;

      // disabled functions
      mcs_lock(const mcs_lock&) = delete;
      mcs_lock(mcs_lock&&) = delete;
      mcs_lock& operator=(const mcs_lock&) = delete;
      mcs_lock& operator=(mcs_lock&&) = delete;
  };

  /// \brief lock using two counters, one for arrival and one for entering time
  /// \details does not implement try_lock
  struct counting_lock
  {
    typedef uint64_t counter_type;

    counting_lock()
    : ctr(0), ticket(0)
    {}

    void lock()
    {
      // hold while other threads are waiting
      while (ticket.val.load(std::memory_order_relaxed) > ctr.val.load(std::memory_order_relaxed) + 1) {}

      counter_type t = ticket.val.fetch_add(1, std::memory_order_relaxed);

      while (t != ctr.val.load(std::memory_order_acquire)) {}
    }

    void unlock()
    {
      counter_type t = ctr.val.load(std::memory_order_relaxed);

      ctr.val.store(t+1, std::memory_order_release);
    }

    bool is_free() const
    {
      return ticket.val.load(std::memory_order_relaxed) == ctr.val.load(std::memory_order_relaxed);
    }

    private:
      aligned_atomic_type<counter_type> ctr;
      aligned_atomic_type<counter_type> ticket;

      // disabled functions
      counting_lock(const counting_lock&) = delete;
      counting_lock(counting_lock&&) = delete;
      counting_lock& operator=(const counting_lock&) = delete;
      counting_lock& operator=(counting_lock&&) = delete;
  };

  //
  // guards and friends

  /// \private
  /// iterate through all tuple elements
  template <size_t I>
  struct tuple_for_all
  {
    template <typename ...Elidable, class F>
    static
    F apply(std::tuple<Elidable&...>& data, F f)
    {
      f(std::get<I-1>(data));
      return tuple_for_all<I-1>::apply(data, f);
    }
  };

  /// \private
  /// base case
  template <>
  struct tuple_for_all<0>
  {
    template <typename ...Elidable, class F>
    static
    F apply(std::tuple<Elidable&...>&, F f)
    {
      return f;
    }
  };

  /// \private
  /// iterate through all tuple elements
  template <size_t I, size_t L>
  struct tuple_find_idx
  {
    template <typename ...Elidable, class F>
    static
    size_t apply(std::tuple<Elidable&...>& data, F f)
    {
      if (f(std::get<I>(data))) return I;
      return tuple_find_idx<I+1, L>::apply(data, f);
    }
  };

  /// \private
  /// base case
  template <size_t I>
  struct tuple_find_idx<I, I>
  {
    template <typename ...Elidable, class F>
    static
    size_t apply(std::tuple<Elidable&...>&, F)
    {
      return I;
    }
  };

  /// \private
  struct is_not_free
  {
    template <class Elidable>
    bool operator()(const Elidable& e)
    {
      return !e.is_free();
    }
  };

  /// \private
  struct acquire_lock
  {
    template <class Lockable>
    void operator()(Lockable& e)
    {
      e.lock();
    }
  };

  /// \private
  struct release_lock
  {
    template <class Lockable>
    void operator()(Lockable& e)
    {
      e.unlock();
    }
  };

#if HTM_ENABLED
  /// \private
  template <class Elideables>
  static inline
  void touch_all(Elideables& elidables)
  {
    static const size_t NUMLOCKS = std::tuple_size<Elideables>::value;

    size_t idx = tuple_find_idx<0, NUMLOCKS>::apply(elidables, is_not_free());

    if (idx < NUMLOCKS) htm::tx::abort<99>(); // abort
  }
#endif

  /// \private
  template <class Lockables>
  static inline
  void lock_all(Lockables& lockables)
  {
    static const size_t NUMLOCKS = std::tuple_size<Lockables>::value;

    tuple_for_all<NUMLOCKS>::apply(lockables, acquire_lock());
  }

  /// \private
  template <class Lockables>
  static inline
  void unlock_all(Lockables& lockables)
  {
    static const size_t NUMLOCKS = std::tuple_size<Lockables>::value;

    tuple_for_all<NUMLOCKS>::apply(lockables, release_lock());
  }

  struct guard_base
  {
    guard_base(const guard_base&)            = delete;
    guard_base& operator=(const guard_base&) = delete;
  };

#if HTM_ENABLED

  /// \private
  template <typename ...Elidable>
  static inline
  bool is_first_free(std::tuple<Elidable&...>& mtx)
  {
    return std::get<0>(mtx).is_free();
  }

  /// \private
  /// creates a lock guard for elidable locks (need to support is_free() )
  ///
  /// \tparam  ...Elidable a sequence of elidable lock types
  /// \details this is an internal class, use the external constructor function below
  template <typename ...Elidable>
  struct elidable_guard : guard_base
  {
      elidable_guard(size_t num_retries, Elidable&... args)
      : elidables(args...)
      {
        while (num_retries > 0 && !htm::tx::begin())
        {
          num_retries = htm::tx::may_retry() ? num_retries-1 : 0;
        }

        if (num_retries)
          touch_all(elidables);
        else
          lock_all(elidables);
      }

      ~elidable_guard()
      {
        if (is_first_free(elidables))
        {
          htm::tx::end();
#if UCL_RUNTIME_DATA
          ++num_commits;
#endif /* UCL_RUNTIME_DATA */
          return;
        }

        unlock_all(elidables);
      }

    private:
      std::tuple<Elidable&...> elidables;
  };

  /// \brief a drop-in replacement for std::lock_guard that elides
  ///        the mutex by default.
  template <class Elidable, size_t NUM_RETRIES = 10>
  struct lock_elision_guard : elidable_guard<Elidable>
  {
    enum { RETRIES = NUM_RETRIES };

    explicit
    lock_elision_guard(Elidable& e)
    : elidable_guard<Elidable>(RETRIES, e)
    {}
  };

  /// creates an elidable lock guard. A guard that attempts to execute a
  ///   a critical section using transactional memory. If the execution does
  ///   not succeed after at most num_tries attempts, it falls back to locked
  ///   execution.
  /// \tparam ...Elidables a sequence of elidable lock types
  ///         elidable lock = mutex + is_free member function
  ///         (e.g., as implemented in lock.hpp)
  /// \param  num_retries number of retries before the locks are taken
  /// \param  a sequence of references to elidable locks
  ///
  /// \details used as follows:
  /// \code
  ///   const auto& guard = ucl::elide_guard(10, m1, m2, m3);
  /// \endcode
  ///
  /// \warning The locks are acquired in sequence. elide_guard DOES NOT
  ///          implement any deadlock prevention.
  template <typename ...Elidables>
  elidable_guard<Elidables...>
  elide_guard(size_t num_retries, Elidables&... args)
  {
    return elidable_guard<Elidables...>(num_retries, args...);
  }

  /// creates an elidable lock-guard, which tries 5 times before falling back
  /// to locks.
  /// \tparam ...Elidables a sequence of elidable lock types
  ///         elidable lock = mutex + is_free member function
  ///         (e.g., as implemented in lock.hpp)
  /// \param  a sequence of references to elidable locks
  ///
  /// \details
  /// \code
  ///   const auto& guard = ucl::elide_guard(m1, m2, m3);
  /// \endcode
  /// \warning The locks are acquired in sequence. elide_guard DOES NOT
  ///          implement any deadlock prevention.
  template <typename ...Elidables>
  elidable_guard<Elidables...>
  elide_guard(Elidables&... args)
  {
    return elidable_guard<Elidables...>(5, args...);
  }
#endif /* HTM_ENABLED */

  /// \private
  /// \brief guards a sequence of locks
  /// \details to avoid naming the Lockable types, use the external
  ///          constructor below.
  template <typename ...Lockable>
  struct lockable_guard : guard_base
  {
      explicit
      lockable_guard(Lockable&... args)
      : lockables(args...)
      {
        lock_all(lockables);
      }

      ~lockable_guard()
      {
        unlock_all(lockables);
      }

    private:
      std::tuple<Lockable&...> lockables;
  };

  /// creates lock-guard for multiple locks at once (similar to std::lock_guard
  ///   but for multiple locks.
  /// \tparam ...M a sequence of mutex types
  /// \param  a sequence of references to locks
  ///
  /// \details used as follows:
  /// \code
  ///   const auto& guard = ucl::lock_guard(m1, m2, m3);
  /// \endcode
  ///
  /// \warning The locks are acquired in sequence. lock_guard DOES NOT
  ///          implement any deadlock prevention.
  template <typename ...M>
  lockable_guard<M...>
  lock_guard(M&... args)
  {
    return lockable_guard<M...>(args...);
  }
}

#endif /* _UNLEASHED_LOCKS_HPP */
