#include <utility>
#include <tbb/task_scheduler_init.h>
#include <tbb/task_group.h>

template <class G>
int task(G& taskgroup)
{
  return 0;
}

int test_tbb()
{
  tbb::task_group g;

  g.run([&g]()->void { task(g); });

  g.wait();
  return 0;
}

int main()
{
  return 0;
}
