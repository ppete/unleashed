
#ifndef _LOCK_SELECTION_HPP
#define _LOCK_SELECTION_HPP 1

#ifdef TEST_MUTEX 
  #include <mutex>
#else
  #include "ucl/spinlock.hpp"
#endif

#if defined TEST_ANDERSON
  #ifndef NUM_THREADS
  #define NUM_THREADS 64
  #endif

  typedef ucl::anderson_lock<NUM_THREADS> default_lock;
#elif defined TEST_TTAS
  typedef ucl::ttas_lock                  default_lock;
#elif defined TEST_TTAS_BO
  typedef ucl::ttas_lock_backoff_default  default_lock;
#elif defined TEST_CLH
  typedef ucl::clh_lock                   default_lock;
#elif defined TEST_MCS
  typedef ucl::mcs_lock                   default_lock;
#elif defined TEST_COUNTING
  typedef ucl::counting_lock              default_lock;
#elif defined TEST_MUTEX
  typedef std::mutex                      default_lock;
#else
  #error "undefined default lock"
#endif /* TEST_* */

#endif /* _LOCK_SELECTION_HPP */
