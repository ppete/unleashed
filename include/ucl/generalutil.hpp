

#ifndef _GENERALUTIL_HPP
#define _GENERALUTIL_HPP

/// avoid unused parameter, variable warnings by calling this function with
///   arguments that are conditionally not used.
template <class... Args>
static inline
void unused(Args&...) {}

#endif
