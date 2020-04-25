

#ifndef _UCL_UNUSED_HPP
#define _UCL_UNUSED_HPP

namespace ucl
{
  /// avoid unused parameter, variable warnings by calling this function with
  ///   arguments that are conditionally not used.
  template <class... Args>
  static inline
  void unused(Args&...) {}
}

#endif /* _UCL_UNUSED_HPP */
