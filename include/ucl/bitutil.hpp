
/// \file bitutil.hpp
/// \brief This file provides a collections of functions to perform
///        bit level manipulations of integers.
///
/// \author Peter Pirkelbauer

#ifndef _BITUTIL_HPP_
#define _BITUTIL_HPP_

#include <cassert>

namespace ucl
{
#ifdef _MSC_VER
  __inline __declspec(naked) int __fastcall bsr32 (size_t)
  { __asm
    { bsr eax, ecx
      ret
    }
  }

#elif defined(__GNUC__)

  /// returns the index of the most significant bit that is set in an unsigned int
  static inline
  int bsr32(unsigned int s)
  {
    assert(s != 0);

    static const int nobits = (sizeof(unsigned int) << 3) - 1;

    // note: use __builtin_clzl, and __builtin_clzll for long and long long
    return nobits - __builtin_clz(s);
  }

#else

#warning "Unsupported Compiler (tested on MVC++, g++, clang++, icc)"

  /// \brief  bit-scan-right: returns the number of the highest bit set in an integer
  /// \param  bnr size_t a non-zero integer value
  /// \pre bnr is not zero
  /// \return the first bit set when bnr is scanned from left to right.
  ///         if bnr == 0, the result is undefined
  static inline
  int bsr32(unsigned int s)
  {
    // \todo optimized implementation to O(log)
    assert(s != 0);

    size_t i = (1 << 31);
    int    cnt = 32;

    while ((i > 0) && ((i & s) == 0))
    {
      i = (i >> 1);
      --cnt;
    }

    return cnt;
  }

#endif /* _MSC_VER, __GNUC__ */


#ifdef __GNUC__

  static inline
  size_t bitcount32(unsigned int n)
  {
    return __builtin_popcount(n);
  }

  static inline
  size_t bitcount64(unsigned long n)
  {
    return __builtin_popcountl(n);
  }

#else /* __GNUC__ */

  /// \cond
  #define BITCOUNT32(x) (((BX_(x)+(BX_(x)>>4)) & 0x0F0F0F0F) % 255)
  #define BX_(x) ((x) - (((x)>>1)&0x77777777) - (((x)>>2)&0x33333333) - (((x)>>3)&0x11111111))
  /// \endcond

  /// \brief  returns the number of bits set in an integer
  /// \param  bnr size_t the integer value
  /// \return the number of bits set in bnr
  static inline
  size_t bitcount32(size_t bnr)
  {
    return BITCOUNT32(bnr);
  }

  /// \cond
  #undef BITCOUNT32
  #undef BX_
  /// \endcond
#endif /* __GNUC__ */
}  // namespace ucl

#endif /* _BITUTIL_HPP_ */
