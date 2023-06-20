#include <cstddef>
#include <omp.h>

int test(int i)
{
  return i;
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
      partialresult = test(i);
    }

    #pragma omp atomic
    total += partialresult;
  }

  return 0;
}

int main()
{
  test_omp(1);
  return 0;
}
