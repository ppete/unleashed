#include <cilk/cilk.h>
#include <cilk/reducer_opadd.h>
#include <cilk/cilk_api.h>

int task(int i)
{
  return 0;
}

size_t partialresult;
#pragma omp threadprivate(partialresult)


int test_omp(int i)
{
  size_t total = 0;

  #pragma omp parallel firstprivate(i) shared(total)
  {
    partialresult = 0;

    #pragma omp single
    #pragma omp taskgroup
    {
      test(i);
    }

    #pragma omp atomic
    total += partialresult;
  }

  return 0;
}

int main()
{
  return 0;
}
