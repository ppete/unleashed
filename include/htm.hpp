/// \brief Hardware Transactional Memory Primitives
///
/// \author  Peter Pirkelbauer
/// \email   pirkelbauer@uab.edu

#ifndef _HTM_HPP

#define _HTM_HPP 1

#include <cassert>

#if defined(__x86_64)
  #include <immintrin.h>
#elif defined(_ARCH_PPC64)
  #include <htmintrin.h>
#else
  #error "Unsupported architecture"
#endif /* defined(arch) */

#ifndef HTM_ENABLED
#define HTM_ENABLED 1
#endif /* HTM_ENABLED */

namespace htm
{
  namespace tx
  {
    /// define transactional states
    enum state { none, active, suspended };

#if defined(__x86_64) // x86

    thread_local size_t txcode = -790614; // some lucky number

    /// starts a transaction
    inline bool begin()
    {
      return ((txcode = _xbegin()) == _XBEGIN_STARTED);
    }

    /// commits a transaction
    inline void end()
    {
      _xend();
    }

    /// aborts a transaction with user code
    /// \tparam X abort code (an unsigned integer value)
    template <unsigned X>
    inline void abort()
    {
      _xabort(X);
    }

    /// returns the abort cause of the last transaction
    inline int abort_cause()
    {
      return txcode;
    }

    /// returns the abort code of the last transaction
    inline int abort_code()
    {
      return _XABORT_CODE(txcode);
    }

    /// returns the transaction's state
    /// \todo check implementations - seems broken when used in elided lock
    inline state transactional()
    {
      if (_xtest()) return none;

      return active;
    }

    /// returns whether a transaction could be retried
    /// note, always true on Intel systems
    inline bool may_retry()
    {
      return true;
    }

    struct intel_x86_tag
    {
      static const size_t len = 80;
    };

    typedef intel_x86_tag arch_tag;
#else // Power8

#ifdef __TM_FENCE__
#define _p8_tbegin(R)  __builtin_tbegin(R);
#define _p8_tend(R)    __builtin_tend(R);
#define _p8_tabort(R)  __builtin_tabort(R);
#else
  /*
   * Enforce compiler barriers on hardware transactions
   * https://sourceware.org/ml/libc-alpha/2016-01/msg00061.html
   *
   * Remove this when glibc drops support for GCC 5.0.
   */
#define _p8_tbegin(R)			\
   ({ __asm__ volatile("" ::: "memory");	\
     unsigned int __ret = __builtin_tbegin (R);	\
     __asm__ volatile("" ::: "memory");		\
     __ret;					\
   })
#define _p8_tabort(R)			\
  ({ __asm__ volatile("" ::: "memory");		\
    unsigned int __ret = __builtin_tabort (R);	\
    __asm__ volatile("" ::: "memory");		\
    __ret;					\
  })
#define _p8_tend(R)			\
   ({ __asm__ volatile("" ::: "memory");	\
     unsigned int __ret = __builtin_tend (R);	\
     __asm__ volatile("" ::: "memory");		\
     __ret;					\
   })
#endif /* __TM_FENCE__  */

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
    inline bool begin()
    {
      return _p8_tbegin(0);
    }

    inline void end()
    {
      _p8_tend(0);
    }

    template <unsigned X>
    inline void abort()
    {
      _p8_tabort(X);
    }
#pragma GCC diagnostic pop

    inline int abort_cause()
    {
      return __builtin_get_texasr();
    }

    inline int abort_code()
    {
      return __builtin_get_texasr();
    }

    inline state transactional()
    {
      unsigned char tx_state = _HTM_STATE (__builtin_ttest ());
      state         res      = none;

      if (tx_state == _HTM_TRANSACTIONAL)
        res = active;
      else if (tx_state == _HTM_SUSPENDED)
        res = suspended;

       return res;
    }

    inline bool may_retry()
    {
      return _TEXASRU_FAILURE_PERSISTENT(__builtin_get_texasru ()) != 0;
    }

    struct power8_tag
    {
      static const size_t len = 21;
    };

    typedef power8_tag arch_tag;
#endif /* x86 / Power8 selection */
  }
} // namespace htm

#ifdef NDEBUG
  #define htm_assert(C) ((void) 0)
#else
  inline
  void htm_assert(bool b)
  {
    if (!b)
    {
      if (htm::tx::transactional() == htm::tx::active)
        htm::tx::end();

      assert(false);
    }
  }
#endif /* NDEBUG */

#endif /* _HTM_HPP */
