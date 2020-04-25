/// \file  pmemory.hpp
/// \brief Memory managers for fine-grained lock based and lock-free data structures
///
/// \details - just_alloc is an allocator that only allocates, but never frees memory.
///          - arena is a memory manager that frees memory at user selected timepoints
///          - pub_scan_manager (PS) is a memory manager that uses
///            publish and scan techniques, similar to M. Michael's hazard pointers
///            and Herlihy et al.'s repeated offender implementation.
///            This technique implies sequentially consistent reads of pointers/elements.
///          - epoch_manager (EP) is a memory manager that uses epochs.
///            EP requires less global synchronization than PS.
///            but does not provide upper bound guarantees on the number of
///            non reclaimable memory. I.e., if a thread terminates in the middle of an
///            operation that uses EP, no memory can be freed.
///          - gc_manager (GC) is a memory manager that works in conjunction with
///            a garbage collector (e.g., Boehm's collector). In essence, GC forwards all
///            calls to the set allocator, which is responsible for releasing the memory
///            when it can no longer be reached. Using GC together with a standard
///            allocator leads to early and HAZARDOUS memory deallocations.
///            To enable the gc_manager, either PMEMORY_GC_ENABLED needs to be defined
///            or Boehm's gc.h must be included before this file.
///
///          The code below specifies the concept for memory managers.
/// \code{.cpp}
/// template <class T, class Alloc>
/// concept memory_manager
/// {
///   // pins an atomic pointer and chooses a mo tag appropriate for the
///   // memory managing technique; The actual memory order is at least as strong.
///   template <std::memory_order mo = std::memory_order_seq_cst>
///   _Tp* pin(std::atomic<_Tp*>& elem);
///
///   // pins an ucl::MarkablePointer and chooses a mo tag appropriate for the
///   //   memory managing technique; The actual memory order is at least as strong.
///   template <std::memory_order mo = std::memory_order_seq_cst, size_t MARKABLE_BITS>
///   typename ucl::MarkablePointer<_Tp, MARKABLE_BITS>::state_type
///   pin(ucl::MarkablePointer<_Tp, MARKABLE_BITS>& elem);
///
///   // pins a non-atomic address (provided to simplify some algorithms)
///   _Tp* pin_addr(_Tp& elem);
///
///   // releases the pointer, and guarantees that all reads/write from/to the
///   //   protected memory have completed.
///   // \details the second argument provides a "hint" to find locate
///   //          the pointer.
///   // \param   mem location whose guard is released
///   // \param   hnt hint, how the location may be guarded
///   void unpin(_Tp* mem, int hnt);
///
///   // releases all pinned pointers, and guarantees that all reads/write
///   //   from/to the protected memory have completed.
///   void unpin_all();
///
///   // begins an operation on shared data
///   // \param sz maximum number of pointers that need to be guarded
///   void beginop(size_t sz);
///
///   // initializes the internal data structure
///   //   if it has not been done yet.
///   // \param sz maximum number of pointers that need to be guarded
///   void initialize_if_needed(size_t sz);
///
///   // ends an operation, and guarantees that all reads/write
///   //   from/to the protected memory have completed.
///   void endop();
///
///   // collects and releases non-referenced memory
///   void release_memory();
///
///   // releases all managed memory unchecked
///   // \pre invoked in quiescent time
///   void qrelease_memory();
///
///   // Releases all memory associated with this memory manager.
///   //   Any thread with access to this memory manager must guarantee
///   //   to only use thread-local resources. Any thread invoking this
///   //   function renders inter thread accesses unusable.
///   // \param mainThread true, if invoked from the main thread   
///   // \pre invoked in quiescent time
///   void qdestruct(bool mainThread);
///
///   // does the memory manager hold unreleased memory?
///   // \details this should be a constant time operation
///   bool has_unreleased_memory() const;
///
///   // returns an estimate of unreleased nodes, if available
///   // \details this may not be a constant time operation
///   size_t count_unreleased_memory() const;
/// };
/// \endcode
///
/// \author
///    Nick Dzugan
///    Amalee Wilson
///    Peter Pirkelbauer

#ifndef _PMEMORY_HPP

#define _PMEMORY_HPP 1



#if defined(GC_H)
  // enable GC if Boehms Collector is included before this header file
  #define PMEMORY_GC_ENABLED 1
#endif

#include <cassert>
#include <atomic>
#include <algorithm>
#include <iterator>
#include <vector>
#include <numeric>
#include <forward_list>
#include <sstream>

#include "ucl/atomicutil.hpp"
#include "ucl/bitutil.hpp"
#include "ucl/unused.hpp"

#ifndef UCL_GC_CXX11_THREAD_CONTEXT
#define UCL_GC_CXX11_THREAD_CONTEXT
namespace {
  /// \brief If GC is enabled and gc-cxx11/gc_cxx11.hpp was included
  ///        gc_cxx_thread_context is defined there.
  ///        Otherwise this introduces a dummy declaration.
  /// \brief dummy thread context when the GC is not used
  struct gc_cxx_thread_context
  {
      /// empty constructor to avoid unused variable warnings
      gc_cxx_thread_context() {}

    private:
      /// should not be allocated on the heap.
      ///   just for documentation, not fool proof.
      void* operator new(size_t size) = delete;
  };
}
#endif /* UCL_GC_CXX11_THREAD_CONTEXT */

/// lockfree code, though not necessarily nonblocking.
namespace lockfree
{
  //
  // memory manager tags

  /// tags memory managers that protect allocation blocks
  struct finegrain {};

  /// tags memory managers that protect an operation
  struct operation {};

  /// tags memory managers that collect memory (e.g., gc_manager and arena)
  struct collected {};

  /// \brief compares allocator tag with a given kind
  /// \tparam T the desired kind
  /// \tparam U the actual kind
  template <class T, class U>
  static inline
  bool is_alloc_kind(U)
  {
    return std::is_same<T, U>::value;
  }

  static const bool FREE_ALWAYS = false; ///< indicates whether the threads should be scanned after
                                         ///  each pointer deallocation.
                                         ///  Of course, this should be false, but interestingly
                                         ///  enough, the overhead seems to be only about 10%
                                         ///  (for 10 threads) when the data is in cache


  /// \private
  /// simple ptr deleter, that uses the provided allocator to release memory.
  template <class _Alloc>
  struct ptr_deleter
  {
    ptr_deleter(_Alloc alloc, typename _Alloc::value_type* elem)
    : node_allocator(alloc), e(elem)
    {}

    ~ptr_deleter()
    {
      node_allocator.deallocate(e, 1);
    }

    _Alloc                             node_allocator;
    typename _Alloc::value_type* const e;
  };


  //
  // Auxiliary classes to manage an allocators PinWall

  /// \brief   Implements an RAII guard to free an allocator's pinwall,
  ///          when a function exits
  /// \tparam  Alloc a lock-free allocator
  template <class Alloc>
  struct guard
  {
      guard(Alloc a, size_t pins)
      : alloc(a)
      {
        alloc.beginop(pins);
      }

      ~guard()
      {
        alloc.endop();
      }

    private:
      Alloc alloc;
  };

  /// \brief   Implements an empty guard for lock-free allocators
  ///          that do not use a pinwall
  struct guardless
  {
    template <class Alloc>
    explicit
    guardless(Alloc, size_t = 0)
    {}
  };

  //
  // auxiliary methods accessing atomic pointers

  /// \private
  template <std::memory_order mo, class _Tp>
  static inline
  _Tp* _ld(std::atomic<_Tp*>& elem)
  {
    return elem.load(mo);
  }

  /// \private
  template <std::memory_order mo, class _Tp, size_t MARKABLE_BITS>
  static inline
  typename ucl::MarkablePointer<_Tp, MARKABLE_BITS>::state_type
  _ld(ucl::MarkablePointer<_Tp, MARKABLE_BITS>& elem)
  {
    //~ return typename ucl::MarkablePointer<_Tp, MARKABLE_BITS>::state_type(nullptr, 0);
    return elem.state(mo);
  }

  /// \private
  template <class _Tp>
  static inline
  _Tp* _ptr(_Tp* ptr)
  {
    return ptr;
  }

  /// \private
  template <class _Tp, class _MARKTYPE>
  static inline
  // _Tp* _ptr(typename ucl::MarkablePointer<_Tp, MARKABLE_BITS>::state_type state)
  _Tp* _ptr(std::pair<_Tp*, _MARKTYPE> state)
  {
    return state.first;
  }

  //
  // lock-free allocator

  /// \private
  /// base class that factors out interfacing with standard allocators
  template <class _Tp, template <class> class _Alloc>
  struct alloc_base
  {
      typedef typename _Alloc<_Tp>::value_type      value_type;
      typedef typename _Alloc<_Tp>::reference       reference;
      typedef typename _Alloc<_Tp>::pointer         pointer;
      typedef typename _Alloc<_Tp>::const_pointer   const_pointer;
      typedef typename _Alloc<_Tp>::const_reference const_reference;
      typedef typename _Alloc<_Tp>::size_type       size_type;
      typedef typename _Alloc<_Tp>::difference_type difference_type;

      // \pp types added for icc
      typedef void*                                 void_pointer;
      typedef const void*                           const_void_pointer;
      typedef std::false_type                       propagate_on_container_copy_assignment;
      typedef std::false_type                       propagate_on_container_move_assignment;
      typedef std::false_type                       propagate_on_container_swap;
      typedef std::false_type                       is_always_equal;

      typedef _Alloc<_Tp>                           base_allocator;
      typedef typename _Alloc<_Tp>::template rebind<void> _AllocVoid;

      explicit
      alloc_base(_Alloc<_Tp> allocator = _Alloc<_Tp>())
      : alloc(allocator)
      {}

      template <class _Up>
      explicit
      alloc_base(const alloc_base<_Up, _Alloc>&)
      : alloc(typename _Alloc<_Up>::template rebind<_Tp>::other())
      {}

      pointer address(reference x) const noexcept
      {
        return alloc.address(x);
      }

      const_pointer address(const_reference x) const noexcept
      {
        return alloc.address(x);
      }

      pointer allocate(size_type n, std::allocator<void>::const_pointer hint = nullptr)
      {
        return alloc.allocate(n, hint);
      }

      template <class U, class... Args>
      void construct (U* p, Args&&... args)
      {
        alloc.construct(p, args...);
      }

      void deallocate(pointer p, size_type n)
      {
        alloc.deallocate(p, n);
      }

      size_type max_size() const noexcept
      {
        return alloc.max_size();
      }

      _Alloc<_Tp> get_allocator()
      {
        return alloc;
      }

    private:
      _Alloc<_Tp>              alloc;
  };

  /// \private
  /// \brief   empty implementation of pin/unpin
  /// \details auxiliary base class for allocators that do not use pin/unpin
  template <class _Tp>
  struct no_pinwall_base
  {
      /// pins an atomic pointer
      template <std::memory_order mo = std::memory_order_seq_cst>
      _Tp* pin(std::atomic<_Tp*>& elem)
      {
        return _ld<mo>(elem);
      }

      /// pins a ucl::MarkablePointer
      template <std::memory_order mo = std::memory_order_seq_cst, size_t MARKABLE_BITS>
      typename ucl::MarkablePointer<_Tp, MARKABLE_BITS>::state_type
      pin(ucl::MarkablePointer<_Tp, MARKABLE_BITS>& elem)
      {
        return _ld<mo>(elem);
      }

      /// pins a non-atomic address (provided to simplify some algorithms)
      _Tp* pin_addr(_Tp& elem)
      {
        return &elem;
      }

      /// releases the pointer
      /// \details the second argument provides a "hint" to find locate
      ///          the pointer.
      void unpin(_Tp*, int) {}

      /// releases all pinned pointers
      void unpin_all() {}

      /// ends an operation
      void endop() {}

      /// collects and releases non-referenced memory
      void release_memory() {}

      /// does not release memory
      void qrelease_memory() {}
      
      /// destructs nothing
      void qdestruct(bool) {}

      /// does the memory manager hold unreleased memory?
      bool has_unreleased_memory() const { return false; }

      /// returns an estimate of unreleased nodes, if available
      size_t count_unreleased_memory() const { return 0; }
  };

  //
  // auxiliary classes

  /// \private
  template <class _Alloc>
  struct Deallocator
  {
    _Alloc alloc;

    explicit
    Deallocator(_Alloc allocator)
    : alloc(allocator)
    {}

    template <class T>
    void operator()(T* ptr)
    {
      alloc.destroy(ptr);

      alloc.deallocate(ptr, 1);
    }
  };

  /// \private
  template <class _Alloc>
  Deallocator<_Alloc> deallocator(_Alloc alloc)
  {
    return Deallocator<_Alloc>{alloc};
  }


  //
  // allocators
  
  /// \brief behaves like a regular allocator that does not manage memory
  /// \details
  ///    frees deallocated memory immediately.
  template <class _Tp, template <class> class _Alloc = std::allocator>
  struct common_alloc : _Alloc<_Tp>
  {
    typedef _Alloc<_Tp> impl;
    
    public:
      // no pinwall
      typedef guardless pinguard;
      typedef collected manager_kind;
      
      using impl::impl;
  };

  /// \brief simple allocator that only allocates BUT DOES NOT FREE
  template <class _Tp, template <class> class _Alloc = std::allocator>
  struct just_alloc : no_pinwall_base<_Tp>, alloc_base<_Tp, _Alloc>
  {
    private:
      typedef alloc_base<_Tp, _Alloc> base;

    public:
      // no pinwall
      typedef guardless  pinguard;
      typedef collected  manager_kind;

      /// \private
      /// rebind to self
      template <class U>
      struct rebind { typedef just_alloc<U, _Alloc> other; };

      /// constructs a just_alloc manager
      explicit
      just_alloc(_Alloc<_Tp> alloc)
      : base(alloc)
      {}

      /// constructs a just_alloc manager
      template <class _Up>
      explicit
      just_alloc(just_alloc<_Up, _Alloc> other)
      : base(other)
      {}

      /// constructs a just_alloc manager
      just_alloc()
      : base()
      {}

      /// disable deallocate
      void deallocate(_Tp*, int) {}
  };


  /// \brief allocator stores pointers to memory until they can be freed
  ///        in a quiescent period.
  template <class _Tp, template <class> class _Alloc = std::allocator>
  struct arena : no_pinwall_base<_Tp>, alloc_base<_Tp, _Alloc>
  {
    private:
      typedef alloc_base<_Tp, _Alloc> base;

      static thread_local std::vector<_Tp> released;

    public:
      // no pinwall
      typedef guardless  pinguard;
      typedef collected  manager_kind;

      /// \private
      /// rebind to self
      template <class U>
      struct rebind { typedef arena<U, _Alloc> other; };

      /// constructs an arena manager
      explicit
      arena(_Alloc<_Tp> alloc)
      : base(alloc)
      {}

      /// constructs an arena manager
      template <class _Up>
      explicit
      arena(arena<_Up, _Alloc> other)
      : base(other)
      {}

      /// constructs an arena manager
      arena()
      : base()
      {}

      /// delay deallocation
      void deallocate(_Tp* el, int)
      {
        released.push_back(el);
      }
      
      /// releases all memory held, under the assumption that the data structure
      ///   is quiescent.
      void qrelease_memory()
      {
        std::for_each(released.begin(), released.end(), deallocator(base::get_allocator()));

        released.clear();
      }

      /// does the memory manager hold unreleased memory?
      bool has_unreleased_memory() const { return released.size() > 0; }

      /// returns an estimate of unreleased nodes, if available
      size_t count_unreleased_memory() const { return released.size(); }
  };

  template <class _Tp, template <class> class _Alloc>
  thread_local
  std::vector<_Tp> arena<_Tp, _Alloc>::released;

#ifdef PMEMORY_GC_ENABLED
  /// \brief memory manager that defers memory reclamation to a
  ///        garbage collector.
  template <class _Tp, template <class> class _Alloc>
  struct gc_manager : no_pinwall_base<_Tp>, alloc_base<_Tp, _Alloc>
  {
    private:
      typedef alloc_base<_Tp, _Alloc> base;

    public:
      // no pinwall
      typedef collected manager_kind;
      typedef guardless pinguard;

      /// \private
      /// rebind to self
      template <class U>
      struct rebind { typedef gc_manager<U, _Alloc> other; };

      /// constructs a gc manager
      explicit
      gc_manager(_Alloc<_Tp> alloc)
      : base(alloc)
      {}

      /// constructs a gc manager
      template <class _Up>
      explicit
      gc_manager(gc_manager<_Up, _Alloc> other)
      : base(other)
      {}

      /// constructs a gc manager
      gc_manager()
      : base()
      {}

      /// disable deallocate and defer to GC
      void deallocate(_Tp*, int) {}
  };

  /// thread context, if GC is enabled
  struct gc_cxx_thread_context
  {
    GC_stack_base sb;

    gc_cxx_thread_context()
    {
      GC_get_stack_base(&sb);
      GC_register_my_thread(&sb);
    }

    ~gc_cxx_thread_context()
    {
      GC_unregister_my_thread();
    }
  };
#endif /* PMEMORY_GC_ENABLED */

  //
  // Auxiliary lock-free stack implementation

  /// \private
  template <class _Tp, class _Alloc>
  struct pub_scan_data
  {
      typedef std::allocator_traits<_Alloc>                          _OrigAlloc_traits;
      typedef typename _OrigAlloc_traits::template rebind_alloc<_Tp> _TpAlloc;
      typedef typename _OrigAlloc_traits::template rebind_alloc<_Tp*> _PpAlloc;

      typedef std::atomic<_Tp*> pinwall_entry;
      typedef std::vector<_Tp*> removal_collection;

      pub_scan_data<_Tp, _Alloc> const * next; ///< next thread's data
      pinwall_entry* const               base; ///< base address of hazard pointers
                                               //  \todo consider to use unique pointer
      std::atomic<pinwall_entry*>        last; ///< memory managed by base
      // size_t                             maxCnt;
      size_t                             cnt;  ///< number of removed ptr since last collection
      size_t                             collection_time; ///< when the next collection should occur
      const size_t                       pinwall_size;
      removal_collection                 rmvd; ///< removed but not yet freed pointers

      // \todo new
      explicit
      pub_scan_data(size_t len)
      : next(nullptr), base(new pinwall_entry[len]), last(base),
        cnt(0), collection_time(len), pinwall_size(len), rmvd()
      {}

      pub_scan_data(const pub_scan_data&) = delete;
      pub_scan_data& operator=(const pub_scan_data&) = delete;

      /// \brief destructor, never to be called as this may generate a hazard
      ///        with a potential concurrent access to the same data
      ~pub_scan_data()
      {
        assert(false); // should never be called
        delete[] base;
      }

      /// \brief internal method that combines handling atomic and markable pointers
      /// \todo choose appropriate memory ordering tags
      template <std::memory_order mo, class _Atomic_Type>
      auto pin_aux(_Atomic_Type& elem) -> decltype(_ld<mo>(elem))
      {
        // Incrementing early, exposes a potentially invalid object
        // This may delay freeing, but should not compromise safety
        pinwall_entry* pwentry = last.load(std::memory_order_relaxed);

        assert(pwentry < base + pinwall_size);
        last.store(pwentry+1, std::memory_order_relaxed); // \mo

        // we confirm later anyhow, so if we do not the most recent value right
        // away we get it later
        auto           lastval = _ld<std::memory_order_relaxed>(elem);
        _Tp*           confirm = nullptr;

        do
        {
          confirm = _ptr(lastval);
          pwentry->store(confirm, std::memory_order_relaxed); // \mo

          // syncs with thread fence in collectAndFree
          std::atomic_thread_fence(std::memory_order_seq_cst); // \mo
          lastval = _ld<std::memory_order_relaxed>(elem);
        } while (confirm != _ptr(lastval)); // \mo

        return lastval;
      }

      /// \brief   pins a non-atomic address
      /// \details allows codes to be more regular.
      ///          e.g., a data structure's entry node is non-modifiable,
      ///                thus does not need to be atomic. Nevertheless, we
      ///                would like to handle the entry node like any other
      ///                atomic node in the data structure.
      _Tp* pin_addr(_Tp& val)
      {
        // Incrementing early, exposes a potentially invalid object
        // This may delay freeing, but should not compromise safety
        pinwall_entry* pwentry = last.load(std::memory_order_relaxed);

        assert(pwentry < base + pinwall_size);
        last.store(pwentry+1, std::memory_order_relaxed);  // \mo

        // no need to confirm, b/c by definition we have no race
        pwentry->store(&val, std::memory_order_relaxed); // \mo
        return &val;
      }

      /// unpins an entry slot range on the pinwall
      void unpin(_Tp const * const elem, int ofs)
      {
        ucl::unused(elem), assert(ofs <= 0);

        // load one past last valid
        //   relaxed sufficient since this is the only writing thread
        pinwall_entry* pwentry = last.load(std::memory_order_relaxed);

        // make ptr point to last valid element
        --pwentry;

        // address of unpinned element
        pinwall_entry* unpinel = pwentry+ofs;

        assert(elem == unpinel->load(std::memory_order_relaxed));

        // move last element into this offset
        if (ofs)
        {
          // sb ordered (single writer)
          _Tp* repl = pwentry->load(std::memory_order_relaxed); // \mo

          unpinel->store(repl, std::memory_order_relaxed); // \mo
        }

        // \todo it seems that we can relax the requirements and just
        //       require last's update to be release instead of a seq_cst fence.
        // \mo   with the update of last the modification to the buffer also
        //       must become visible.
        std::atomic_thread_fence(std::memory_order_seq_cst); // \mo

        last.store(pwentry, std::memory_order_relaxed); // \mo release?
        assert(pwentry >= base);
      }

      /// unpins an entry slot range on the pinwall
      void check(_Tp const * const elem)
      {
        pinwall_entry* const zz = last.load(std::memory_order_relaxed);

        ucl::unused(elem), assert(std::find(base, zz, elem) != zz);
      }

      void unpin_all()
      {
        last.store(base, std::memory_order_relaxed);
      }

      void endop()
      {
        unpin_all();
      }

      bool needs_collect() const { return collection_time <= cnt; }
  };

  /// \private
  /// \brief Iterator that goes through a hazard pointer data stack
  template <class _ScanData>
  struct scan_iterator : std::iterator<std::forward_iterator_tag, const _ScanData>
  {
    typedef std::iterator<std::forward_iterator_tag, const _ScanData> base;

    explicit
    scan_iterator(_ScanData* where)
    : pos(where)
    {}

    bool operator==(const scan_iterator& that) const
    {
      return this->pos == that.pos;
    }

    bool operator!=(const scan_iterator& that) const
    {
      return this->pos != that.pos;
    }

    scan_iterator operator=(const scan_iterator& that)
    {
      this->pos = that.pos;
      return *this;
    }

    const typename base::value_type&
    operator*()
    {
      return *pos;
    }

    scan_iterator& operator++()
    {
      assert(pos != nullptr);

      pos = pos->next;
      return *this;
    }

    scan_iterator operator++(int)
    {
      scan_iterator tmp(*this);

      assert(pos != nullptr);
      pos = pos->next;

      return tmp;
    }

    const _ScanData* pos;
  };

  /// \private
  /// \brief Hazard Pointer Collector
  template <class T, class _Alloc>
  struct pub_scan_scanner
  {
    typedef typename pub_scan_data<T, _Alloc>::removal_collection removal_collection;

    removal_collection collectedPointers;
    const size_t       maxlen;

    pub_scan_scanner(size_t sz, size_t maxPtrs)
    : collectedPointers(), maxlen(sz)
    {
      collectedPointers.reserve(maxPtrs * maxlen);
    }

    void operator()(const pub_scan_data<T, _Alloc>& data)
    {
      // \pp acquire needed to sync with rel in unpin
      const std::atomic<T*>* limit = data.last.load(std::memory_order_acquire); // \mo
      const std::atomic<T*>* base  = data.base;

      // to remain consistent with unpin operations (which unpin a slot by moving the
      // pointer from the last position into the slot of the unpinned pointer),
      // we traverse from back to front.
      // Suppose that we traverse from the base to limit:
      // The collector thread could be at position 3, when an unpin operation
      // at location 2 moves the last element into 2 and then remove the last
      // element. This way, the collector would never see 2.
      // By traversing the list in reverse order:
      //   (a) we have either seen the last ptr before it was moved into 2
      //   (b) we have not seen it, in which case the ptr was added after we
      //       started scanning and can be ignored
      while (limit != base)
      {
        --limit;
        // T* ptr = limit->load(std::memory_order_relaxed); // \mo
        T* ptr = limit->load(); // \mo

        collectedPointers.push_back(ptr);
      }
    }

    removal_collection result()
    {
      return std::move(collectedPointers);
    }
  };

  /// \private
  /// \brief convenience function to create Hazard Pointer Collectors
  template <class T, class _Alloc>
  pub_scan_scanner<T, _Alloc> pubScanScanner(size_t sz, size_t maxPtrs)
  {
    return pub_scan_scanner<T, _Alloc>(sz, maxPtrs);
  }

  /// \private
  /// \brief   functor testing whether a ptr is NOT in a less-than sorted
  ///          sequence of pointers
  /// \details the ptr pointers also need to be supplied in a less-than
  ///          sorted form
  /// \return  true, if ptr is not in cont; false otherwise
  template <class _Tp, class _Alloc>
  struct Pinned
  {
    typedef std::vector<_Tp*, _Alloc> container;

    typename container::const_iterator const aa;
    typename container::const_iterator const zz;

    explicit
    Pinned(const container& cont)
    : aa(cont.begin()), zz(cont.end())
    {
    }

    bool operator()(_Tp const * const ptr)
    {
      return std::binary_search(aa, zz, ptr);
    }
  };

  /// \private
  template <template <class, class> class _Cont, class _Tp, class _Alloc>
  Pinned<_Tp, _Alloc> pinned(_Cont<_Tp*, _Alloc> const & cont)
  {
    return Pinned<_Tp, _Alloc>(cont);
  }

  /// \private
  static inline
  size_t threshold(size_t len)
  {
    return len * ucl::bsr32(len+1);
  }

  /// \brief   Allocator implementing a publish and scan memory manager
  ///          (aka Hazard pointers).
  /// \warning This publish and scan implementation is not reentrant
  /// \note    The current implementation is unsuitable for programming styles
  ///          that frequently create new threads.
  template <class _Tp, template <class> class _Alloc = std::allocator>
  struct pub_scan_manager : alloc_base<_Tp, _Alloc>
  {
    private:
      typedef alloc_base<_Tp, _Alloc> base;

      // choose the RECLAMATION_PERIOD based on maximum number
      //   of threads and pinwall size.
      //   N = |threads| * pinwallsize
      //   RECLAMATION_PERIOD > N lg N
      // static const size_t RECLAMATION_PERIOD = 2048;

    public:
      using base::construct;
      using base::max_size;
      using base::address;

      typedef finegrain                 manager_kind;
      typedef typename base::value_type value_type;

      /// \private
      template <class U>
      struct rebind { typedef pub_scan_manager<U, _Alloc> other; };

      typedef std::atomic<value_type*>                  pinwall_type;
      typedef pinwall_type*                             pinwall_base;
      typedef std::atomic<pinwall_base>                 pinwall_pointer;
      typedef guard<pub_scan_manager<_Tp, _Alloc> >     pinguard;

      /// constructs a publish and scan manager
      explicit
      pub_scan_manager(_Alloc<_Tp> alloc)
      : base(alloc)
      {}

      template <class _Up>
      explicit
      pub_scan_manager(pub_scan_manager<_Up, _Alloc> other)
      : base(other)
      {}

      /// constructs a publish and scan manager
      pub_scan_manager()
      : base()
      {}

    private:
      /// collects all currently published pointers from available threads
      typename pub_scan_data<value_type, _Alloc<_Tp> >::removal_collection
      collectAllPointersFromThreads()
      {
        typedef pub_scan_data<value_type, _Alloc<_Tp> >    pub_scan_data;
        typedef typename pub_scan_data::removal_collection removal_collection;
        typedef pub_scan_scanner<value_type, _Alloc<_Tp> > pub_scanner;
        typedef scan_iterator<pub_scan_data>               pub_scan_iterator;

        pub_scan_data*     curr = allPinWalls.load(std::memory_order_relaxed);
        pub_scan_iterator  pos(curr);
        pub_scan_iterator  end(nullptr);
        removal_collection res = std::for_each(pos, end, pub_scanner(pinWall->collection_time, (size_t)1)).result();

        return res;
      }

    public:
      /// \brief  runs the scan and collect algorithm and frees memory that can no
      ///         longer be referenced by other threads.
      void release_memory()
      {
        typedef pub_scan_data<value_type, _Alloc<_Tp> >    pub_scan_data;
        typedef typename pub_scan_data::removal_collection removal_collection;

        assert(pinWall);

        // syncs with thread fence in pin
        std::atomic_thread_fence(std::memory_order_seq_cst);

        // collect and sort all currently pinned pointers
        removal_collection                    pinnedPtrs = collectAllPointersFromThreads();

        std::sort(pinnedPtrs.begin(), pinnedPtrs.end());

        // sort previously freed pointers
        removal_collection&                   rmvd = pinWall->rmvd;
        typename removal_collection::iterator aa   = rmvd.begin();
        typename removal_collection::iterator zz   = rmvd.end();

        // std::sort(aa, zz);

        // make the non pinned pointers available to be freed
        typename removal_collection::iterator pos = std::partition(aa, zz, pinned(pinnedPtrs));

        // free memory that is no longer used
        std::for_each(pos, zz, deallocator(base::get_allocator()));

        // eliminate already freed from set
        rmvd.erase(pos, zz);

        // reset iteration count
        pinWall->cnt = 0;
        pinWall->collection_time = pinWall->pinwall_size + threshold(pinnedPtrs.size());
      }

      /// releases all memory under assumption that data structure is quiescent
      /// i.e., no other thread has access to these nodes
      void qrelease_memory()
      {
        typedef pub_scan_data<value_type, _Alloc<_Tp> >    pub_scan_data;
        typedef typename pub_scan_data::removal_collection removal_collection;

        removal_collection&                   rmvd = pinWall->rmvd;

        std::for_each(rmvd.begin(), rmvd.end(), deallocator(base::get_allocator()));
        rmvd.clear();
      }
      
      void qdestruct(bool mainThread)
      {
        delete pinWall;
        pinWall = nullptr;
        
        if (mainThread) allPinWalls.store(nullptr, std::memory_order_relaxed);
      }

      /// returns true if this thread has unreleased memory
      bool has_unreleased_memory() const
      {
        return (count_unreleased_memory() > 0);
      }

      /// returns the number of unreleased memory blocks
      size_t count_unreleased_memory() const
      {
        assert(pinWall);

        return pinWall->rmvd.size();
      }

      /// allocates a new memory block
      value_type*
      allocate(size_t num, const void* hint=0)
      {
        assert(num == 1);

        value_type* res = base::allocate(num, hint);
        return res;
      }

      /// destroy and deallocate the object pointed by obj
      ///   Memory management delays destruction
      ///   until no other thread holds a reference to the same object
      void deallocate(value_type* obj, size_t num)
      {
        assert(pinWall);
        ucl::unused(num), assert(num == 1); // currently the allocator is limited to a single object

        pinWall->rmvd.push_back(obj);
        ++pinWall->cnt;

        if (FREE_ALWAYS || pinWall->needs_collect()) release_memory();
      }

      void initialize_if_needed(size_t szPinWall)
      {
        typedef pub_scan_data<value_type, _Alloc<_Tp> > pub_scan_data;

        if (!pinWall)
        {
          // \todo new
          pinWall = new pub_scan_data(szPinWall);

          pub_scan_data* curr = allPinWalls.load(std::memory_order_relaxed);

          do
          {
            pinWall->next = curr;
          } while (!allPinWalls.compare_exchange_weak(curr, pinWall, std::memory_order_release, std::memory_order_relaxed)); // \mo
        }
      }

      /// \brief starts an operation and creates slots for szPinWall entries
      void beginop(size_t szPinWall)
      {
        initialize_if_needed(szPinWall);
      }

      /// \brief pins an atomic pointer
      template <std::memory_order mo = std::memory_order_seq_cst>
      value_type*
      pin(std::atomic<value_type*>& elem)
      {
        assert(pinWall);
        return pinWall->template pin_aux<mo>(elem);
      }

      /// \brief pins a ucl::MarkablePointer
      template <std::memory_order mo = std::memory_order_seq_cst, size_t MARKABLE_BITS>
      typename ucl::MarkablePointer<value_type, MARKABLE_BITS>::state_type
      pin(ucl::MarkablePointer<value_type, MARKABLE_BITS>& elem)
      {
        assert(pinWall);
        return pinWall->template pin_aux<mo>(elem);
      }

      /// \brief pins a non-atomic address
      // \todo  choose appropriate memory ordering tags
      // \todo  consider renaming / removing
      value_type*
      pin_addr(value_type& val)
      {
        assert(pinWall);
        return pinWall->pin_addr(val);
      }

      /// unpins an entry slot range on the pinwall
      void unpin(value_type const * const elem, int ofs)
      {
        assert(pinWall);
        pinWall->unpin(elem, ofs);
      }

      /// unpins all entries
      void unpin_all()
      {
        assert(pinWall);
        pinWall->unpin_all();
      }

      /// ends an operation
      void endop()
      {
        assert(pinWall);
        pinWall->endop();
      }

    private:
      /// creates new storage for a Thread, unless the
      typedef std::atomic<pub_scan_data<value_type, _Alloc<_Tp> >*> pinwall_entry;

      static pinwall_entry                                          allPinWalls;
      static thread_local pub_scan_data<value_type, _Alloc<_Tp> >*  pinWall;
  };

  /// \private
  template <class _Tp, template <class> class _Alloc>
  typename pub_scan_manager<_Tp, _Alloc>::pinwall_entry
  pub_scan_manager<_Tp, _Alloc>::allPinWalls(nullptr);

  /// \private
  template <class _Tp, template <class> class _Alloc>
  thread_local
  pub_scan_data<typename pub_scan_manager<_Tp, _Alloc>::value_type, _Alloc<_Tp> >*
  pub_scan_manager<_Tp, _Alloc>::pinWall(nullptr);

  //
  // epoch management scheme

  /// \private
  /// \brief stores pointers collected at a given epoch
  template <class _T, class _Alloc>
  struct release_entry
  {
    // rebind alloc to release entry types
    typedef typename _Alloc::template rebind<size_t>::other size_t_alloc;
    typedef typename _Alloc::template rebind<_T*>::other    tp_alloc;
    typedef std::vector<size_t>               epoch_vector;
    typedef std::vector<_T*>                      epoch_ptr_vector;

    epoch_vector     epochtime;
    epoch_ptr_vector epochptrs;

    explicit
    release_entry(size_t timesz = 0)
    : epochtime(), epochptrs()
    {
      epochtime.reserve(timesz);
      epochptrs.reserve(timesz);
    }
  };

  /// \private
  /// \brief stores all collected epochs until they can be freed
  template <class _Tp, class _Alloc>
  struct epoch_data
  {
    typedef std::forward_list<release_entry<_Tp, _Alloc> > removal_collection;

    epoch_data<_Tp, _Alloc> const * next;  ///< next thread's data
    std::atomic<size_t>             epoch; ///< odd numbers indicate busy epochs

    ucl::aligned_type<size_t>       dummy;
#if UCL_RUNTIME_DATA
    size_t                          cnt_collects = 0;
#endif /* UCL_RUNTIME_DATA */
    removal_collection              rmvc;


    bool needs_collection()
    {
      assert(!rmvc.empty());

      typename removal_collection::value_type& entry = rmvc.front();

      return entry.epochtime.size() + 1 < entry.epochptrs.size();
    }

    epoch_data()
    : next(nullptr), epoch(0), dummy(), rmvc(1)
    {
      assert(!rmvc.empty());
    }
  };

  /// \private
  /// test whether an epoch counter indicates activity
  static inline
  bool active(size_t epoch)
  {
    return epoch & 1;
  }


  /// \private
  /// \brief finds the epochs that can be reclaimed
  template <class _Tp, class _Alloc>
  struct passed_epoch_finder
  {
    typedef epoch_data<_Tp, _Alloc>                   epoch_data_t;
    typedef release_entry<_Tp, _Alloc>                release_entry_t;
    typedef typename epoch_data_t::removal_collection removal_collection;
    typedef typename release_entry_t::epoch_vector    epoch_vector;

    // true if freeepoch and currepoch produce a "mismatch".
    //   mismatch =
    static bool epoch_passed(size_t freeepoch, size_t currepoch)
    {
      assert(freeepoch <= currepoch);

      return (  freeepoch != currepoch // unchanged epoch component
             || !active(freeepoch)     // active
             );
    }

    typename epoch_vector::const_reverse_iterator const curraa;

    explicit
    passed_epoch_finder(const epoch_vector& vec)
    : curraa(vec.rbegin())
    {}

    bool operator()(typename removal_collection::value_type& entry)
    {
      typename epoch_vector::reverse_iterator const entryzz = entry.epochtime.rend();

      // mismatch finds the first component where the current epoch (curraa) has
      //   not advanced passed the epoch when elements where freed
      //   (entry.epochtime).
      //   - reverse iterator, b/c epochs are recorded in reverse order, with
      //     the oldest component being at the last position.
      //   - the length of the current epoch vector may be longer than the epoch
      //     vector at the time when objects were freed, thus we iterator over
      //     entry.epochtime.
      typename epoch_vector::reverse_iterator pos = std::mismatch(entry.epochtime.rbegin(), entryzz, curraa, epoch_passed).first;

      return entryzz == pos;
    }
  };


  /// \private
  template <class _FwdIter, class _UnaryPred>
  _FwdIter fwd_find(_FwdIter before_aa, _FwdIter zz, _UnaryPred pred)
  {
    assert( before_aa != zz);

    _FwdIter pos = before_aa;

    ++pos;
    while ( pos != zz && !pred(*pos) )
    {
      ++pos; ++before_aa;
    }

    return before_aa;
  }


  /// \private
  template <class _Alloc>
  struct epoch_deallocator
  {
     typedef typename std::allocator_traits<_Alloc>::value_type value_type;

     _Alloc alloc;

     explicit
     epoch_deallocator(_Alloc allocator)
     : alloc(allocator)
     {}

     void operator()(release_entry< value_type, _Alloc >& epoch)
     {
       std::for_each( epoch.epochptrs.begin(),
                      epoch.epochptrs.end(),
                      [this](value_type* ptr)
                      {
                        alloc.destroy(ptr);
                        alloc.deallocate(ptr,1);
                      }
                    );
     }
  };

  /// \private
  template <class _Alloc>
  epoch_deallocator<_Alloc>
  epochDeallocator(_Alloc alloc)
  {
    return epoch_deallocator<_Alloc>(alloc);
  }

  /// \brief   epoch based memory manager (compare RCU methods)
  /// \details Each thread keeps its own epoch counter. The epoch counter is
  ///          incremented when an operation starts and again incremented when
  ///          an operation ends.
  ///          Memory can be reclaimed, when a thread is non active (epoch
  ///          counter is even) or a thread's epoch counter was advanced since
  ///          the memory was made inaccessible.
  ///          Hence, the epoch_manager does not provide an upper bound on
  ///          non-reclaimable memory blocks. I.e., if a thread fails bad while
  ///          an operation is active - w/o exception handling being involved -,
  ///          no more memory can be reclaimed.
  /// \warning The current implementation is not reentrant.
  /// \note    The current implementation is unsuitable for programming styles
  ///          that frequently create new threads.
  template <class _Tp, template <class> class _Alloc = std::allocator >
  struct epoch_manager : alloc_base<_Tp, _Alloc>
  {
      typedef operation                   manager_kind;
      typedef alloc_base<_Tp, _Alloc>     base;
      typedef typename base::value_type   value_type;
      typedef typename base::pointer      pointer;
      typedef typename base::size_type    size_type;

      /// no pinwall
      typedef guard<epoch_manager<value_type, _Alloc> > pinguard;

      /// \private
      /// rebind to self
      template <class U>
      struct rebind { typedef epoch_manager<U, _Alloc> other; };

      using base::construct;
      using base::max_size;
      using base::address;

      /// constructs a publish and scan manager with the supplied allocator
      explicit
      epoch_manager(_Alloc<_Tp> alloc)
      : base(alloc)
      {}

      template <class _Up>
      epoch_manager(const epoch_manager<_Up, _Alloc>& orig)
      : base(orig)
      {}

      /// constructs a publish and scan manager
      epoch_manager()
      : base()
      {}

      pointer allocate(size_type n, std::allocator<void>::const_pointer hint = nullptr)
      {
        return base::allocate(n, hint);
      }

      /// safely loads the content of an atomic variable
      template<std::memory_order mo = std::memory_order_seq_cst>
      value_type* pin(std::atomic<value_type*>& elem)
      {
        return _ld<mo>(elem);
      }

      /// safely loads the content of ucl::MarkablePointer
      template<std::memory_order mo = std::memory_order_seq_cst, size_t MARKBITS>
      typename ucl::MarkablePointer<value_type, MARKBITS>::state_type
      pin(ucl::MarkablePointer<value_type, MARKBITS>& elem)
      {
        return _ld<mo>(elem);
      }

      /// just returns a pointer to the element
      /// \details available to handle normal pointers and shared pointers uniformly
      value_type* pin_addr(value_type& elem)
      {
        return &elem;
      }

      /// unpin not needed
      void unpin(value_type*, int) {}

      void
      collect_epochs(typename release_entry<value_type, _Alloc<_Tp> >::epoch_vector& res)
      {
        typedef scan_iterator<epoch_data<value_type, _Alloc<_Tp> > > epoch_iterator;
        typedef epoch_data<value_type, _Alloc<_Tp> >                 epoch_data_t;

        epoch_data_t*   curr = allepochData.load(std::memory_order_relaxed);
        epoch_iterator  end(nullptr);

        // load all epochs
        for (epoch_iterator pos(curr); pos != end; ++pos)
        {
          // \mo we load relaxed; the global barrier must have been seen before
          res.push_back((*pos).epoch.load(std::memory_order_relaxed));
        }
      }

      void
      push_curr_list(typename epoch_data<value_type, _Alloc<_Tp> >::removal_collection::value_type&& epochdesc)
      {
        epochdata->rmvc.front().swap(epochdesc);
      }

      /// \private
      void free_unreferenced_memory()
      {
        typedef epoch_data<value_type, _Alloc<_Tp> >          epoch_data_t;
        typedef typename epoch_data_t::removal_collection     removal_collection;
        typedef typename removal_collection::iterator         removal_iterator;
        typedef release_entry<value_type, _Alloc<_Tp> >       release_entry_t;
        typedef typename release_entry_t::epoch_vector        epoch_vector;
        typedef passed_epoch_finder<value_type, _Alloc<_Tp> > passed_epoch_finder_t;

        removal_collection&    rmvc   = epochdata->rmvc;
        epoch_vector&          epoch  = rmvc.front().epochtime;
        removal_iterator const preaa  = rmvc.before_begin();
        removal_iterator const zz     = rmvc.end();

        // could use binary search (if we had a vector), but it is expected
        //   that we do not keep too many removal records around.
        removal_iterator const prepos = fwd_find(preaa, zz, passed_epoch_finder_t(epoch));
        removal_iterator       pos    = prepos;

        std::advance(pos, 1);

        std::for_each(pos, zz, epochDeallocator(base::get_allocator()));
        rmvc.erase_after(prepos, zz);
      }

      /// \brief  compares epochs of releasable elements with current epochs of threads.
      ///         Frees all memory that can no longer be referenced by other threads.
      void release_memory()
      {
        typedef release_entry<value_type, _Alloc<_Tp> > release_entry_t;

        assert(epochdata);

#if UCL_RUNTIME_DATA
        ++epochdata->cnt_collects;
#endif /* UCL_RUNTIME_DATA */

        // syncs with thread fence in the epoch counter
        std::atomic_thread_fence(std::memory_order_seq_cst);

        // collect all epochs and add current epoch to list
        release_entry_t& curr_entry = epochdata->rmvc.front();

        collect_epochs(curr_entry.epochtime);

        // get current clock size
        const size_t     current_clock_size = curr_entry.epochtime.size();

        // use the epochs to clear old memory
        free_unreferenced_memory();

        // add a new entry for the next epoch
        epochdata->rmvc.emplace_front( current_clock_size );
      }

      /// \brief  compares epochs of releasable elements with current epochs of threads.
      ///         Frees all memory that can no longer be referenced by other threads.
      void release_memory_if_needed()
      {
        if (epochdata->needs_collection())
        {
          release_memory();
          // epochdata->collepoch += 8*epochdata->reclaimation_increment();
          // std::cout << epochdata->collepoch << std::endl;
        }
      }

      /// releases all memory
      void qrelease_memory()
      {
        assert(epochdata);

        typedef epoch_data<value_type, _Alloc<_Tp> >      epoch_data_t;
        typedef typename epoch_data_t::removal_collection removal_collection;

        removal_collection& rmvc = epochdata->rmvc;

        std::for_each(rmvc.begin(), rmvc.end(), epochDeallocator(base::get_allocator()));

        rmvc.clear();
        rmvc.emplace_front( 0 );

        assert(!rmvc.empty());
      }
      
      void qdestruct(bool mainThread)
      {
        delete epochdata;
        epochdata = nullptr;
        
        if (mainThread) allepochData.store(nullptr, std::memory_order_relaxed); 
      }

      /// \brief  returns whether releasable memory blocks were held back
      ///         due to other threads being still able to reference them.
      bool has_unreleased_memory() const
      {
        typedef epoch_data<value_type, _Alloc<_Tp> >      epoch_data_t;
        typedef typename epoch_data_t::removal_collection removal_collection;

        typename removal_collection::iterator aa = epochdata->rmvc.begin();
        typename removal_collection::iterator zz = epochdata->rmvc.end();

        return (  (aa != zz && (++aa) != zz)
               || (!epochdata->rmvc.front().epochptrs.empty())
               );
      }

      /// \brief  returns the number of memory blocks that are held back due to other threads
      ///         being still able to reference them.
      size_t count_unreleased_memory() const
      {
         return std::accumulate( epochdata->rmvc.begin(),
                                 epochdata->rmvc.end(),
                                 0,
                                 [](size_t sz, release_entry<value_type, _Alloc<_Tp> >& entry)
                                 {
                                   return sz + entry.epochptrs.size();
                                 }
                               );
      }

      /// unpins all pointers (a noop for epoch manager)
      void unpin_all() {}

      /// increments the epoch counter
      void endop()
      {
        size_t epochval = epochdata->epoch.load(std::memory_order_relaxed) + 1;

        assert(!active(epochval)); // exiting an operation

        // \mo release in order to prevent reordering with preceding operations.
        epochdata->epoch.store(epochval, std::memory_order_release);

        release_memory_if_needed();
        assert(!epochdata->rmvc.empty());
      }

      void initialize_if_needed(size_t = 0)
      {
        typedef epoch_data<value_type, _Alloc<_Tp> > epoch_data_t;

        if (!epochdata)
        {
          // \todo \new
          epochdata = new epoch_data_t();

          epoch_data_t* curr = allepochData.load(std::memory_order_relaxed);

          do
          {
            epochdata->next = curr;
          } while (!allepochData.compare_exchange_weak(curr, epochdata, std::memory_order_release, std::memory_order_relaxed)); // \mo
          return ;
        }
      }

      /// \brief starts an operation and creates a data entry for the epoch managers if needed
      /// \param dummy parameter to provide consistent interface with other managers
      void beginop(size_t)
      {
        initialize_if_needed();

        // \mo this thread updates its own epoch, thus we can use relaxed
        size_t epochval = epochdata->epoch.load(std::memory_order_relaxed) + 1;

        assert(active(epochval)); // entering an operation
        epochdata->epoch.store(epochval, std::memory_order_seq_cst);
        assert(!epochdata->rmvc.empty());
      }

      void deallocate(value_type* entry, int)
      {
        assert(epochdata);
        assert(!epochdata->rmvc.empty());

        epochdata->rmvc.front().epochptrs.push_back(entry);
      }

    // private:
      /// creates new storage for a Thread, unless the
      typedef std::atomic<epoch_data<value_type, _Alloc<_Tp> >*> epoch_wall_type;
      typedef epoch_data<value_type, _Alloc<_Tp> >               epoch_data_rec_t;

      static epoch_wall_type                                     allepochData;
      static thread_local epoch_data<value_type, _Alloc<_Tp>>*   epochdata;
  };

  template <class _Tp, template <class> class _Alloc>
  typename epoch_manager<_Tp, _Alloc>::epoch_wall_type
  epoch_manager<_Tp, _Alloc>::allepochData(nullptr);

  template <class _Tp, template <class> class _Alloc>
  thread_local
  epoch_data<typename epoch_manager<_Tp, _Alloc>::value_type, _Alloc<_Tp> >*
  epoch_manager<_Tp, _Alloc>::epochdata(nullptr);

#if UCL_RUNTIME_DATA
  template <class Alloc>
  static inline
  void log_epoch_telemetry(std::ostream& os, Alloc)
  {
    typedef typename Alloc::epoch_data_rec_t epoch_data_rec_t;

    const epoch_data_rec_t* curr = Alloc::allepochData.load();
    size_t        x              = 0;
    size_t        total_collects = 0;

    while (curr != nullptr)
    {
      total_collects += curr->cnt_collects;

      os << ": e" << (++x)
         << "  [" << curr->cnt_collects << "]"
         << std::endl;

      curr = curr->next;
    }

    os << ":::: etotals: " << total_collects
       << std::endl;
  }
#endif /* UCL_RUNTIME_DATA */

}

#endif /* _PMEMORY_HPP */
