
// tests compilability of MarkablePointer

#include <iostream>
#include <atomic>
#include "atomicutil.hpp"

int main()
{
  typedef ucl::MarkablePointer<int>::mark_type mark_t;

  ucl::MarkablePointer<int> mp(nullptr);

  const int* ptr1 = mp.load();
  const int* ptr2 = mp.load(std::memory_order_relaxed);

  assert(ptr1 == ptr2 && ptr1 == nullptr);

  int   base = 1;

  mp.store(&base, 1);

  std::cout << "mp = " << *mp << std::endl;

  int     other = 2;
  int*    currptr = &base;
  mark_t  currmark = 0;

  const bool res1 = mp.compare_exchange_strong(currptr, &other, currmark, 1);
  std::cout << "mp = " << *mp << "  " << res1 << std::endl;

  currmark = 1;
  const bool res2 = mp.compare_exchange_weak(currptr, &other, currmark, 0);
  std::cout << "mp = " << *mp << "  " << res2 << std::endl;
  std::cout << "marks = " << mp.mark() << std::endl;

  bool res3 = mp.mark_weak(&other, 1);
  std::cout << "marks = " << mp.mark() << " " << res3 << std::endl;

  bool res4 = mp.mark_strong(&other, 0);
  std::cout << "ptr = " << mp.state().first << "  mark = " << mp.state().second << " " << res4 << std::endl;

  return 0;
}
