#ifndef _UCL_THREAD
#define _UCL_THREAD 1

#include <thread>

namespace ucl
{
#if defined(__cpp_lib_jthread)
  using thread = std::jthread;
#else
  using thread = std::thread;
#endif /* __has_cpp_attribute */

}

#endif /* _UCL_THREAD */
