#include <cilk/cilk.h>
#include <cilk/reducer_opadd.h>
#include <cilk/cilk_api.h>

int task()
{
  return 0;
}

int test_cilk()
{
  cilk_spawn task();
  return 0;
}

int main()
{
  return 0;
}
