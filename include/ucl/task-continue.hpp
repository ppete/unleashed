#ifndef _UCL_TASK_CONTINUE_HPP
#define _UCL_TASK_CONTINUE_HPP 1

#include <cassert>
#include <atomic>
#include <functional>

#include "ucl/unused.hpp"

namespace ucl
{
  //
  // model continuations

  /// \private
  /// \brief class holding counter + continuation function
  struct count_fun_base
  {
    typedef int_fast32_t counter_t;

    explicit
    count_fun_base(counter_t c)
    : cnt(c)
    {}
    
    virtual ~count_fun_base() = default;
    virtual void execute()    = 0;

    std::atomic< counter_t > cnt;
  };
  
  template <class Fn>
  struct count_fun : count_fun_base
  {
    using count_fun_base::counter_t;

    count_fun(Fn&& fn, counter_t c)
    : count_fun_base(c), fun(std::move(fn))
    {}

    ~count_fun() override {}
    
    void execute() override { fun(); }

    Fn                       fun;
    std::atomic< counter_t > cnt;
  };

  /// \brief   models a task continuation
  /// \details counts the number of references to the continuation.
  ///          The last active reference executes the lambda function as
  ///          continuation. Continuations can be copied, thus the number
  ///          of references can change during execution.
  ///          - Example 1: empty continuation
  /// \code
  ///   continuation nil;
  /// \endcode
  ///
  ///          - Example 2: merge sort using continuations. outer is the parent's
  ///            continuation, which gets executed, after the inner continuation
  ///            completes.
  /// \code
  ///   continuation inner([low, mid, high, outer]()->{ merge(low,mid,high); });
  /// \endcode
  struct continuation
  {
    typedef count_fun_base::counter_t counter_t;

    public:
      /// ctor to implement "manual move"
      /// \details
      ///   used in benchmarks to implement non-waiting OpenMP tasks
      explicit
      continuation(count_fun_base* existing)
      : cont(existing)
      {
        //~ std::cerr << "ctor(*)" << std::endl;
      }

      /// constructs a new continuation with an initial count
      template <class Fn>
      explicit
      continuation(Fn fun, counter_t initial_cnt)
      : cont(new count_fun<Fn>(std::move(fun), initial_cnt))
      {
        assert(initial_cnt > 0);
        //~ std::cerr << "ctor" << std::endl;
      }

      /// constructs an empty continuation
      continuation()
      : cont(nullptr)
      {
        //~ std::cerr << "ctor()" << std::endl;
      }

      /// moves a continuation without changing the count
      continuation(continuation&& other)
      : cont(other.cont)
      {
        //~ std::cerr << "mctor" << std::endl;
        other.cont = nullptr;
      }

      /// copies a continuation and updates the count accordingly
      continuation(const continuation& other)
      : cont(other.cont)
      {
        //~ std::cerr << "cctor" << std::endl;
        if (cont) inc();
      }

      /// move-assigns a continuation without changing the count
      continuation& operator=(continuation&& other)
      {
        //~ std::cerr << "m=" << std::endl;
        continuation tmp(std::move(*this));

        cont = other.cont;
        other.cont = nullptr;
        return *this;
      }

      /// copies a continuation and increases the count
      continuation& operator=(const continuation& other)
      {
        //~ std::cerr << "c=" << std::endl;
        continuation tmp(std::move(*this));

        if (other.cont)
        {
          cont = other.cont;
          inc();
        }

        return *this;
      }

      /// copies the continuation without changing the count
      ///   \note if several references are created from the same
      ///         origin, the change in references can be tracked
      ///         externally and updated only once.
      ///         However, care needs to be taken that the count does
      ///         not reach 0 preliminary, otherwise the continuation
      ///         gets executed and released too early.
      ///         The examples/tasks/health/health.cc benchmark contains
      ///         an example.
      continuation copy_external_count() const
      {
        return continuation(cont);
      }

      /// releases num references that were externally counted
      void release_external(counter_t num)
      {
        assert(cont != nullptr);

        const counter_t numrefs = dec(num);

        if (numrefs == num)
        {
          std::atomic_thread_fence(std::memory_order_acquire);
          cont->execute();
          delete cont;
        }

        cont = nullptr;
      }

      /// releases the continuation, if no external references were created.
      /// \param num must match the number of outstanding references
      void release_local(counter_t num)
      {
        assert(cont);
        assert(cont->cnt.load(std::memory_order_relaxed) == num);
        ucl::unused(num);

        cont->execute();
        delete cont;
        cont = nullptr;
      }

      /// releases a reference to the continuation
      ///   if this was the last reference, the cont. will get executed.
      ~continuation()
      {
        //~ std::cerr << "dtor" << std::endl;
        if (cont != nullptr) release_external(1);
      }

      count_fun_base* omp_ptr() const
      {
        return cont;
      }

#if 0
      /// \private
      /// \note experimental
      count_fun_base* unlink()
      {
        count_fun* res = cont;

        cont = nullptr;
        return res;
      }
#endif
    private:
      void inc()
      {
        assert(cont);

        // \mo no need to acquire anything besides task
        cont->cnt.fetch_add(1, std::memory_order_relaxed);
      }

      size_t dec(counter_t x = 1)
      {
        assert(cont);

        // last one does neither need to decrement nor release
        if (cont->cnt.load(std::memory_order_relaxed) == x)
        {
          return x;
        }

        // \mo potentially release data computed by task
        return cont->cnt.fetch_sub(x, std::memory_order_release);
      }

      count_fun_base* cont;

    //
    // static elements

    public:
      static constexpr
      counter_t counter_max()
      {
        // return max / 2 to leave room for both increments and decrements
        return std::numeric_limits<counter_t>::max() / 2;
      }
  };


#if 0
  template <class T>
  struct shared_data {};

  /// \brief   models a shared pointer w/ some optimizations for tasks built-in
  template <class T>
  struct shared_ptr
  {
    typedef typename shared_data<T>::counter_t counter_t;

    public:
      /// ctor to implement "manual move"
      /// \details
      ///   used in benchmarks to implement non-weaiting OpenMP tasks
      explicit
      shared_ptr(shared_data<T>* existing)
      : mem(existing)
      {
        //~ std::cerr << "ctor(*)" << std::endl;
      }

      /// constructs a new continuation with an initial count
      explicit
      shared_ptr(T* ptr, counter_t initial_cnt = 1)
      : mem(new shared_data<T>(ptr, initial_cnt))
      {
        assert(initial_cnt > 0);
        assert(cont->cnt.load(std::memory_order_relaxed) == initial_cnt);
        //~ std::cerr << "ctor" << std::endl;
      }

      /// constructs an empty ptr
      shared_ptr()
      : mem(nullptr)
      {
        //~ std::cerr << "ctor()" << std::endl;
      }

      /// moves a pointer without changing the count
      shared_ptr(shared_ptr&& other)
      : mem(other.mem)
      {
        //~ std::cerr << "mctor" << std::endl;
        other.mem = nullptr;
      }

      /// move-assigns a continuation without changing the count
      shared_ptr& operator=(shared_ptr&& other)
      {
        //~ std::cerr << "m=" << std::endl;
        continuation tmp(std::move(*this));

        cont = other.cont;
        other.cont = nullptr;
        return *this;
      }

      /// copies a pointer and updates the count accordingly
      shared_ptr(const shared_ptr& other) = delete;

      /// copies a continuation and increases the count
      shared_ptr& operator=(const shared_ptr& other) = delete;

      /// copies the continuation without changing the count
      ///   \note if several references are created from the same
      ///         origin, the change in references can be tracked
      ///         externally and updated only one.
      ///         However, care needs to be taken that the count does
      ///         not reach 0 preliminarily, otherwise the continuation
      ///         gets executed and released too early.
      shared_ptr copy_external_count() const
      {
        return shared_ptr(mem);
      }

      /// releases num references that were externally counted
      void release_external(counter_t num)
      {
        assert(cont != nullptr);

        const counter_t numrefs = dec(num);

        if (numrefs == num)
        {
          std::atomic_thread_fence(std::memory_order_acquire);
          delete data;
        }

        cont = nullptr;
      }

      /// releases the continuation, if no external references were created.
      /// \param num must match the number of outstanding references
      void release_local(counter_t num)
      {
        assert(cont);
        assert(cont->cnt.load(std::memory_order_relaxed) == num);
        ucl::unused(num);

        cont->fun();
        delete data;
        cont = nullptr;
      }


      /// releases a reference to the continuation
      ///   if this was the last reference, the cont. will get executed.
      ~continuation()
      {
        //~ std::cerr << "dtor" << std::endl;
        if (data != nullptr) release_external(1);
      }

      count_fun* omp_ptr() const
      {
        return cont;
      }

    private:
      void inc()
      {
        assert(false);
        assert(cont);
        //~ ++numatomics;
        // \mo no need to acquire anything besides task
        cont->cnt.fetch_add(1, std::memory_order_relaxed);
      }

      size_t dec(counter_t x = 1)
      {
        assert(cont);

        if (cont->cnt.load(std::memory_order_release) == x)
          return x;

        //~ ++numatomics;
        // \mo potentially release data computed by task
        return cont->cnt.fetch_sub(x, std::memory_order_release);
      }

      shared_mem<T*>* data;

    //
    // static elements

    public:
      static constexpr
      counter_t counter_max()
      {
        return std::numeric_limits<counter_t>::max();
      }

      //~ static size_t numatomics;
  };
#endif /* NOT YET */
} // end namespace ucl

#endif /* _UCL_TASK_CONTINUE_HPP */
