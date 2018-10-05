
#ifndef COMMON_QTHREADS_HPP
#define COMMON_QTHREADS_HPP 1

#include <sstream>

template <class V>
static inline
void qthr_setenv(std::string name, V val)
{
  std::stringstream str;

  str << name << "=" << val;

  size_t numchr = str.str().size()+1;
  char*  envset = static_cast<char*>(malloc(sizeof(char) * numchr));

  memcpy(envset, str.str().c_str(), numchr);
  putenv(envset);

  // putenv may use the string, thus it MUST NOT be freed
  // free envset;
}

static inline
void init_qthreads(size_t numthreads, size_t stacklen = 0)
{
  assert(numthreads);

  size_t worker_per_shepherd = 2;
  size_t num_shepherds = (numthreads+worker_per_shepherd-1) / worker_per_shepherd;

  qthr_setenv("QTHREAD_INFO", 2); // print config output
  qthr_setenv("QTHREAD_HWPAR", numthreads);
  qthr_setenv("QTHREAD_NUM_SHEPHERDS", num_shepherds);

  if (stacklen) qthr_setenv("QTHREAD_STACK_SIZE", stacklen);

  qthread_initialize();
  // qthread_init(num_shepherds);
}

#endif /* COMMON_QTHREADS_HPP */
