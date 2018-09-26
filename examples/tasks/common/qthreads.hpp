
#ifndef COMMON_QTHREADS_HPP
#define COMMON_QTHREADS_HPP 1

#include <sstream>

template <class V>
static inline
void qthr_setenv(std::string name, V val)
{
  std::stringstream str;

  str << name << "=" << val;

  char* envset = new char[str.str().size()+1];

  memcpy(envset, str.str().c_str(), str.str().size()+1);
  putenv(envset);

  // delete[] envset;
}

static inline
void init_qthreads(size_t numthreads, size_t stacklen = 0)
{
  assert(numthreads);

  size_t num_shepherds = (numthreads+1) / 2;

  qthr_setenv("QTHREAD_INFO", 2); // print config output
  qthr_setenv("QTHREAD_HWPAR", numthreads);
  qthr_setenv("QTHREAD_NUM_SHEPHERDS", num_shepherds);

  if (stacklen) qthr_setenv("QTHREAD_STACK_SIZE", stacklen);

  qthread_initialize();
  // qthread_init(num_shepherds);
}

#endif /* COMMON_QTHREADS_HPP */
