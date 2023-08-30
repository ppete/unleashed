/// \file htm.hpp
/// \brief Abstracts hardware transactional memory primitives for x86, Power, and ARM
/// \author  Peter Pirkelbauer
#ifndef _HTM_HPP
#define _HTM_HPP 1

#include <cassert>

#if defined(__x86_64)
  #include <immintrin.h>
#elif defined(_ARCH_PPC64)
  #include <htmintrin.h>
#elif defined(__aarch64__)
  #include <arm_acle.h>
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
    /// \todo check implementations 
    inline state status()
    {
      if (_xtest()) return none;

      return active;
    }

    /// returns whether a transaction could be retried
    /// note, always true on Intel systems
    inline constexpr
    bool may_retry()
    {
      return true;
    }

    struct intel_x86_tag
    {
      enum : size_t { len = 80 };
    };

    typedef intel_x86_tag arch_tag;
#elif defined (_ARCH_PPC64) // Power

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

    inline state status()
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
      enum : size_t { len = 21 };
    };

    typedef power8_tag arch_tag;
#else /* ARM */
    thread_local uint64_t txstate = 8273827;

    inline bool begin()
		{
			return ((txstate = __tstart()) == 0);
		}

    inline void end()
    {
		  __tcommit();
		}

		template <unsigned X>
		inline void abort()
		{
			__tcancel(X);
		}

		inline state status()
		{
	    return __test() ? none : active;  
    }

		inline bool may_retry()
		{
			return (txstate & _TMFAILURE_RTRY) == _TMFAILURE_RTRY; 
		}

		struct arm_tag
		{
			enum : size_t { len = 21 };
		};

		typedef arm_tag arch_tag;

#endif /* x86 / Power8 selection / ARM */
  }
} // namespace htm

#ifdef NDEBUG
  #define htm_assert(C) ((void) 0)
#else
  static inline
  void htm_assert(bool b)
  {
    if (!b)
    {
      if (htm::tx::status() != htm::tx::none)
        htm::tx::end();

      assert(false);
    }
  }
#endif /* NDEBUG */

#endif /* _HTM_HPP */
