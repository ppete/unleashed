/**
 * A simple reduction object for task-frameworks that do not support reductions
 * naturally.
 *
 * NOTE, THE DESIGN IS ONLY SUITABLE FOR BENCHMARKING, BUT WOULD BREAK EASILY
 * IN REAL APPLICATIONS!!!
 */

#ifndef _SIMPLE_REDUCER_HPP
#define _SIMPLE_REDUCER_HPP

#include <cassert>

#include "atomicutil.hpp"

namespace uab
{
  /// DO NOT USE IN REAL CODE!!!
  template <class T, size_t MAXTHR = 256>
  struct simple_reducer
  {
    typedef uab::aligned_type<T, CACHELINESZ> entry_t;

    simple_reducer(T init = T())
    : ctr(1)
    {
      assert(thrid == MAXTHR);

      for (size_t i = 0; i < MAXTHR; ++i)
      {
        data[i].val = init;
      }

      thrid = 0;
    }

    T& get_object()
    {
      if (thrid == MAXTHR)
      {
        thrid = ctr.fetch_add(1, std::memory_order_relaxed);
        assert(thrid < MAXTHR);
      }

      return data[thrid].val;
    }

    T get_value()
    {
      for (size_t i = 1; i < MAXTHR; ++i)
      {
        data[0].val += data[i].val;
      }

      return std::move(data[0].val);
    }

    simple_reducer& operator+=(const T& rhs)
    {
      get_object() += rhs;
      return *this;
    }

    entry_t             data[MAXTHR];
    std::atomic<size_t> ctr;

    static thread_local size_t thrid; /// MULTIPLE REDUCERS FAIL HERE!
  };

  template <class T, size_t MAXTHR>
  thread_local size_t simple_reducer<T,MAXTHR>::thrid = MAXTHR;
}

#endif /* _SIMPLE_REDUCER_HPP */
