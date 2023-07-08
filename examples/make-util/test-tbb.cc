#include <utility>

#if __has_include(<tbb/task_scheduler_init.h>) && __has_include(<tbb/task_group.h>)

#include <tbb/task_scheduler_init.h>
#include <tbb/task_group.h>

namespace TBB = tbb;

#elif __has_include(<oneapi/tbb/global_control.h>) && __has_include(<oneapi/tbb/task_group.h>)

#include <oneapi/tbb/global_control.h>
#include <oneapi/tbb/task_group.h>

//~ oneapi::tbb::global_control global_limit(oneapi::tbb::global_control::max_allowed_parallelism, 2);

namespace TBB = oneapi::tbb;

#else

#error "TBB header files not found."

#endif

template <class G>
int task(G& taskgroup)
{
  return 0;
}

void test_tbb()
{
  TBB::task_group g;

  g.run([&g]()->void { task(g); });
  g.wait();
}

int main()
{
  test_tbb();
  return 0;
}
