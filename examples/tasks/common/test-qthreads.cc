
#include <cstdlib>
#include <sstream>
#include <qthread/qthread.hpp>
#include <qthread/sinc.h>

aligned_t task(void*)
{
  return 0;
}

int test_qthreads()
{
  aligned_t ret;

  qthread_fork(task, nullptr, &ret);
  return 0;
}

void init_qthreads()
{
  std::stringstream str;

  str << "QTHREAD_HWPAR=" << 4;

  char* envset = new char[str.str().size()+1];

  memcpy(envset, str.str().c_str(), str.str().size()+1);
  putenv(envset);

  qthread_initialize();
}

int main()
{
  init_qthreads();
  test_qthreads();

  return 0;
}
