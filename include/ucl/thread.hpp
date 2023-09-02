#ifndef _UCL_THREAD
#define _UCL_THREAD 1

#include <thread>

namespace ucl
{
#if __cplusplus >= 202002L
  using thread = std::jthread;
#else /* !C++20 */
  using thread = std::thread;
#endif /* !C++20 */
}

#endif /* _UCL_THREAD */
