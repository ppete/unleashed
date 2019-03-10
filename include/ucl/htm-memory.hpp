/// \file htm-memory.hpp
/// \brief Memory managers for data structures relying on HTM
///
/// \details - just_alloc is an allocator that only allocates, but never frees memory.
///          - pub_scan_manager (PS) is a memory manager that uses
///            publish and scan techniques, similar to M. Michael's hazard pointers
///            and Herlihy et al.'s repeated offender implementation.
///          - ref_counter is reference counting scheme for HTM
///          - epoch_manager is a coarse grain scheme that HTM makes non-blocking
/// \author  Peter Pirkelbauer

#ifndef _HTMMEMORY_HPP
#define _HTMMEMORY_HPP 1

#if defined(GC_H)
  // enable GC if Boehms Collector is included
  #define HTM_MEMORY_GC_ENABLED true
#endif

#include <cassert>
#include <atomic>
#include <algorithm>
#include <iterator>
#include <vector>
#include <numeric>
#include <forward_list>
#include <sstream>

#include <iostream>

#include "atomicutil.hpp"
#include "bitutil.hpp"
#include "htm.hpp"


#ifndef UCL_GC_CXX11_THREAD_CONTEXT
#define UCL_GC_CXX11_THREAD_CONTEXT
/// \brief If GC is enabled and gc-cxx11/gc_cxx11.hpp was included
///        gc_cxx_thread_context is defined there.
///        Otherwise this introduces a dummy declaration.
namespace
{
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

namespace htm
{
  struct finegrain {};
  struct operation {};
  struct collected {};

  template <template <typename> class _Alloc, class _Tp>
  struct rebinder
  {
    typedef typename _Alloc<_Tp>::template rebind<_Tp> type;
  };

  static const bool FREE_ALWAYS = 0; ///< indicates whether the threads should be scanned after
                                     ///  each pointer deallocation.
                                     ///  Of course, this should be false, but interestingly
                                     ///  enough, the overhead seems to be only about 10%
                                     ///  (for 10 threads) when the data is in cache


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

      bool active() { return _active(&alloc); }

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
    guardless(Alloc, size_t)
    {}

    bool active() { return true; }
  };

  //
  // auxiliary methods accessing atomic pointers

  template <std::memory_order mo, class _Tp>
  static inline
  _Tp* _ld(std::atomic<_Tp*>& elem)
  {
    return elem.load(mo);
  }

  template <std::memory_order mo, class _Tp, size_t MARKABLE_BITS>
  static inline
  typename ucl::MarkablePointer<_Tp, MARKABLE_BITS>::state_type
  _ld(ucl::MarkablePointer<_Tp, MARKABLE_BITS>& elem)
  {
    //~ return typename ucl::MarkablePointer<_Tp, MARKABLE_BITS>::state_type(nullptr, 0);
    return elem.state(mo);
  }

  template <class _Tp>
  static inline
  _Tp* _ptr(_Tp* ptr)
  {
    return ptr;
  }

  template <class _Tp, class _MARKTYPE>
  static inline
  // _Tp* _ptr(typename ucl::MarkablePointer<_Tp, MARKABLE_BITS>::state_type state)
  _Tp* _ptr(std::pair<_Tp*, _MARKTYPE> state)
  {
    return state.first;
  }

  //
  // lock-free allocator

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

      // pointer allocate(size_type n, typename _AllocVoid::const_pointer hint = nullptr)
      pointer allocate(size_type n, std::allocator<void>::const_pointer hint = nullptr)
      {
        return alloc.allocate(n, hint);
      }

      template <class U, class... Args>
      void construct (U* p, Args&&... args)
      {
        alloc.construct(p, args...);
      }

      void deallocate (pointer p, size_type n)
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

  /// \brief   empty implementation of pin/unpin
  /// \details auxiliary base class for allocators that do not use pin/unpin
  template <class _Tp>
  struct no_pinwall_base
  {
      /// pins a range of pointers
      template <class _TrivialIterator>
      void txpin(_TrivialIterator, _TrivialIterator, void*) {}

      /// releases all pinned pointers
      void unpin_all() {}

      /// ends an operation
      void endop() {}

      /// collects and releases non-referenced memory
      void release_memory() {}

      size_t has_unreleased_memory() { return 0; }
  };

  struct stats_data
  {
    stats_data* next;
    size_t      stats;
    size_t      statsctr;
    size_t      collectctr;
  };

  /// \brief Iterator that goes through a thread data stack
  template <class _ScanData>
  struct scan_iterator : std::iterator<std::forward_iterator_tag, _ScanData>
  {
    typedef std::iterator<std::forward_iterator_tag, _ScanData> base;

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

    typename base::value_type&
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

    _ScanData* pos;
  };



  /// \brief simple allocator that only allocates BUT DOES NOT FREE
  template <class _Tp, template <class> class _Alloc = std::allocator>
  struct just_alloc : no_pinwall_base<_Tp>, alloc_base<_Tp, _Alloc>
  {
    private:
      typedef alloc_base<_Tp, _Alloc> base;

    public:
      // no pinwall
      // typedef guardless  pinguard;
      typedef guard<just_alloc<_Tp, _Alloc> > pinguard;
      typedef collected  manager_kind;

      // rebind to self
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

      ///
      void beginop(size_t)
      {
        if (!justdata)
        {
          justdata = new stats_data;

          stats_data* curr = alljusts.load(std::memory_order_relaxed);

          do
          {
            justdata->next = curr;
          } while (!alljusts.compare_exchange_weak(curr, justdata, std::memory_order_release, std::memory_order_relaxed)); // \mo
        }
      }


      scan_iterator<stats_data>
      statsbegin() { return scan_iterator<stats_data>(alljusts); }

      scan_iterator<stats_data>
      statslimit() { return scan_iterator<stats_data>(nullptr); }

      void stats(size_t num)
      {
        justdata->stats += num;
        ++justdata->statsctr;
      }

      static thread_local stats_data*  justdata;
      static std::atomic<stats_data*>  alljusts;
  };

  template <class _Tp, template <class> class _Alloc>
  std::atomic<stats_data*>
  just_alloc<_Tp, _Alloc>::alljusts(nullptr);

  template <class _Tp, template <class> class _Alloc>
  thread_local
  stats_data*
  just_alloc<_Tp, _Alloc>::justdata(nullptr);


  //
  // Auxiliary classes

  template <class _Tp, class _Alloc>
  struct alignas(CACHELINESZ) pub_scan_data
  {
      typedef std::allocator_traits<_Alloc>                           _OrigAlloc_traits;
      typedef typename _OrigAlloc_traits::template rebind_alloc<_Tp*> _TpAlloc;

      typedef std::atomic<_Tp*>           pinwall_entry;
      typedef std::vector<_Tp*, _TpAlloc> removal_collection;

      pub_scan_data<_Tp, _Alloc>* next;            ///< next thread's data
      pinwall_entry* const        base;            ///< base address of hazard pointers
                                                   //  \todo consider to use unique pointer
      std::atomic<pinwall_entry*> last;            ///< memory managed by base
      size_t                      delocctr;        ///< number of removed ptr since last collection
      size_t                      collection_time; ///< when the next collection should occur
      const size_t                pinwall_size;
      size_t                      stats;
      size_t                      statsctr;
      size_t                      collectctr;
      removal_collection          rmvd;            ///< removed but not yet freed pointers

      explicit
      pub_scan_data(size_t len)
      : next(nullptr), base(new pinwall_entry[len]), last(base),
        delocctr(0), collection_time(len), pinwall_size(len),
        stats(0), statsctr(0), collectctr(0), rmvd()
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

      /// pins a range of pointers
      template <class _ForwardIterator>
      void txpin(_ForwardIterator aa, _ForwardIterator zz)
      {
        pinwall_entry* pwentry = base;

        while (aa != zz)
        {
          htm_assert(pwentry < base + pinwall_size);
          pwentry->store(*aa, std::memory_order_relaxed);
          ++pwentry;
          ++aa;
        }

        last.store(pwentry, std::memory_order_relaxed);
      }

      void unpin_all()
      {
        last.store(base, std::memory_order_relaxed);
      }

      void endop()
      {
        unpin_all();
      }
  };


  template <class T>
  static inline
  bool needs_collect(const T t)
  {
    return t->collection_time <= t->delocctr;
  }


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

  /// \brief convenience function to create Hazard Pointer Collectors
  template <class T, class _Alloc>
  pub_scan_scanner<T, _Alloc> pubScanScanner(size_t sz, size_t maxPtrs)
  {
    return pub_scan_scanner<T, _Alloc>(sz, maxPtrs);
  }

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

  template <template <class, class> class _Cont, class _Tp, class _Alloc>
  Pinned<_Tp, _Alloc> pinned(_Cont<_Tp*, _Alloc> const & cont)
  {
    return Pinned<_Tp, _Alloc>(cont);
  }


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

  template <class _Alloc>
  Deallocator<_Alloc> deallocator(_Alloc alloc)
  {
    return Deallocator<_Alloc>{alloc};
  }

  template <class T>
  static inline
  void param_unused(const T&) {}

  static inline
  size_t threshold(size_t len)
  {
    return len * ucl::bsr32(len+1);
  }

  /// \brief   Allocator implementing a publish and scan memory manager
  /// \warning This publish and scan implementation is not reentrant
  /// \note    The current implementation is unsuitable for programming styles
  ///          that frequently create new threads
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

        // make the non pinned pointers available to be freed
        typename removal_collection::iterator pos = std::partition(aa, zz, pinned(pinnedPtrs));

        // free memory that is no longer used
        std::for_each(pos, zz, deallocator(base::get_allocator()));

        // eliminate already freed from set
        rmvd.erase(pos, zz);

        // reset iteration count
        pinWall->delocctr = 0;
        pinWall->collection_time = 2 * pinWall->pinwall_size + threshold(pinnedPtrs.size());
        ++pinWall->collectctr;
      }

      /// returns the number of unreleased memory blocks
      size_t has_unreleased_memory()
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
        assert(num == 1), param_unused(num); // currently the allocator is limited to a single object

        pinWall->rmvd.push_back(obj);
        ++pinWall->delocctr;

        if (FREE_ALWAYS || needs_collect(pinWall)) release_memory();
      }

      /// \brief starts an operation and creates slots for szPinWall entries
      void beginop(size_t szPinWall)
      {
        typedef pub_scan_data<value_type, _Alloc<_Tp> > pub_scan_data;

        if (!pinWall)
        {
          pinWall = new pub_scan_data(szPinWall);

          pub_scan_data* curr = allPinWalls.load(std::memory_order_relaxed);

          do
          {
            pinWall->next = curr;
          } while (!allPinWalls.compare_exchange_weak(curr, pinWall, std::memory_order_release, std::memory_order_relaxed)); // \mo
        }
      }

      template <class _ForwardIterator>
      void txpin(_ForwardIterator aa, _ForwardIterator zz, void*)
      {
        htm_assert(pinWall);
        pinWall->txpin(aa, zz);
      }

      /// unpins all entries
      void unpin_all()
      {
        htm_assert(pinWall);
        pinWall->unpin_all();
      }

      /// ends an operation
      void endop()
      {
        assert(pinWall);
        pinWall->endop();
      }

      void stats(size_t num)
      {
        pinWall->stats += num;
        ++pinWall->statsctr;
      }

      scan_iterator< pub_scan_data<value_type, _Alloc<_Tp> > >
      statsbegin()
      {
        return scan_iterator< pub_scan_data<value_type, _Alloc<_Tp> > >(allPinWalls);
      }

      scan_iterator< pub_scan_data<value_type, _Alloc<_Tp> > >
      statslimit()
      {
        return scan_iterator< pub_scan_data<value_type, _Alloc<_Tp> > >(nullptr);
      }

      //

    private:
      /// creates new storage for a Thread, unless the
      typedef std::atomic<pub_scan_data<value_type, _Alloc<_Tp> >*> pinwall_entry;

      static pinwall_entry                                          allPinWalls;
      static thread_local pub_scan_data<value_type, _Alloc<_Tp> >*  pinWall;
  };

  template <class _Tp, template <class> class _Alloc>
  typename pub_scan_manager<_Tp, _Alloc>::pinwall_entry
  pub_scan_manager<_Tp, _Alloc>::allPinWalls(nullptr);

  template <class _Tp, template <class> class _Alloc>
  thread_local
  pub_scan_data<typename pub_scan_manager<_Tp, _Alloc>::value_type, _Alloc<_Tp> >*
  pub_scan_manager<_Tp, _Alloc>::pinWall(nullptr);


  //
  // Reference Counting

  typedef stats_data ref_data;

  template <class _Tp>
  struct counted : _Tp
  {
    template <class... Args>
    explicit
    counted(Args&&... args)
    : _Tp(args...), dummy(-1), open(0), closed(0)
    {}

    alignas(CACHELINESZ) size_t dummy;
    alignas(CACHELINESZ) size_t open;
    alignas(CACHELINESZ) size_t closed;
  };

  template <class _Tp>
  static
  bool non_removable(counted<_Tp>* el)
  {
    // std::cout << el->open << "  " << el->closed << std::endl;
    return el->open != el->closed;
  }

  /// \brief   Allocator implementing a reference counting scheme
  /// \warning This publish and scan implementation is not reentrant
  /// \note    The current implementation is unsuitable for programming styles
  ///          that frequently create new threads
  template <class _Tp, template <class> class _Alloc = std::allocator>
  struct ref_counter : alloc_base<counted<_Tp>, _Alloc>
  {
    private:
      typedef alloc_base<counted<_Tp>, _Alloc>                       base;

    public:
      using base::construct;
      using base::address;

      typedef finegrain                 manager_kind;
      typedef typename base::value_type counted_type; ///< internal value_type
      typedef counted_type value_type;
      // typedef _Tp                       value_type;   ///< external value_type
                                                      // \todo ^^^^ is this in line w/ the standard?

      template <class U>
      struct rebind { typedef ref_counter<U, _Alloc> other; };

      // conversely to dynamic collect algorithms, the pinwall is not read
      //    by other threads, thus we do not need atomic<>.
      typedef counted_type*                      pinwall_type;
      typedef pinwall_type*                      pinwall_base;
      typedef guard<ref_counter<_Tp, _Alloc> >   pinguard;

    private:
      typedef std::allocator_traits<_Alloc<counted<_Tp> > >                    _OrigAlloc_traits;

    public:
      typedef typename _OrigAlloc_traits::template rebind_alloc<counted_type*> _PtrAlloc;

    private:
      static thread_local pinwall_base                          last;
      static thread_local pinwall_base                          pinw;
      static thread_local std::vector<counted_type*, _PtrAlloc> rmvd;
      static thread_local ref_data*                             refdata;
      static std::atomic<ref_data*>                             allrefs;

    public:
      /// constructs a publish and scan manager
      explicit
      ref_counter(_Alloc<counted<_Tp> > alloc)
      : base(alloc)
      {}

      template <class _Up>
      explicit
      ref_counter(ref_counter<_Up, _Alloc> other)
      : base(other)
      {}

      /// constructs a publish and scan manager
      ref_counter()
      : base()
      {}

      /// pins a range of pointers [ ptr[aa]; ptr[zz] )
      ///   and unpins all no longer needed ptrs.
      template <class _BidirectionalIterator>
      void txpin(_BidirectionalIterator aa, _BidirectionalIterator zz, pinwall_type pinned)
      {
        htm_assert(pinw);

        pinwall_type* l = pinw;

        // we run reverse to overlap with the entries on pinwall
        //   \todo regularize the implementation
        --aa; --zz;

        // a) unpin all objs that are no longer needed
        // b) pin new objs
        // c) skip already pinned objs
        while ((l != last) && (aa != zz))
        {
          if (*zz == *l)
          {
            // aa is already pinned through current entry
            pinned = (*l);
            ++l;
          }
          else if (*zz != pinned)
          {
            // aa is not pinned
            //  unpin l, and pin aa
            pinned = static_cast<counted_type*>(*zz);
            ++(pinned->open);
            ++((*l)->closed);
            (*l) = pinned;
            ++l;
          }
          else
          {
            // aa is already pinned through previous entry -> no update
            // l is not overwritten -> no update
          }

          --zz;
        }

        // unpin remaining old objs
        pinwall_type* l1 = l;

        while (l1 != last)
        {
          ++((*l1)->closed);

          ++l1;
        }

        // pin remaining new pointers
        while (aa != zz)
        {
          if (*zz != pinned)
          {
            pinned = (*l) = static_cast<counted_type*>(*zz);
            ++(pinned->open);
            ++l;
          }

          --zz;
        }

        // set new last element
        last = l;
      }

      /// releases all pinned pointers
      void unpin_all()
      {
        static const size_t len = tx::arch_tag::len;
        size_t              i   = len;

        while (true)
        {
          if (tx::begin())
          {
            while (last != pinw && i > 0)
            {
              --i;
              --last;
              ++((*last)->closed);
            }

            tx::end();
            if (last == pinw) return;
          }
        }
      }

      /// ends an operation
      void endop() { unpin_all(); }

      /// collects and releases non-referenced memory
      void release_memory()
      {
        auto aa      = rmvd.begin();
        auto zz      = rmvd.end();

        // make the non pinned pointers available to be freed
        auto pos     = std::partition(aa, zz, non_removable<_Tp>);

        // free memory that is no longer used
        std::for_each(pos, zz, deallocator(base::get_allocator()));

        // eliminate already freed from set
        rmvd.erase(pos, zz);

        ++refdata->collectctr;
      }

      size_t has_unreleased_memory() { return rmvd.size() > 0; }

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
        assert(pinw && num == 1), unused(num); // currently the allocator is limited to a single object

        rmvd.push_back(static_cast<counted_type*>(obj));
        if (rmvd.size() >= 1024) release_memory();
      }

      /// \brief starts an operation and creates slots for szPinWall entries
      void beginop(size_t szPinWall)
      {
        typedef counted_type* counted_type_ptr;

        if (!pinw)
        {
          last = pinw = new counted_type_ptr[szPinWall];
          refdata = new ref_data;

          ref_data* curr = allrefs.load(std::memory_order_relaxed);

          do
          {
            refdata->next = curr;
          } while (!allrefs.compare_exchange_weak(curr, refdata, std::memory_order_release, std::memory_order_relaxed));
        }
      }

      void stats(size_t num)
      {
        refdata->stats += num;
        ++refdata->statsctr;
      }

      scan_iterator< ref_data >
      statsbegin()
      {
        return scan_iterator< ref_data >(allrefs);
      }

      scan_iterator< ref_data >
      statslimit()
      {
        return scan_iterator< ref_data >(nullptr);
      }
  };


  template <class _Tp, template <class> class _Alloc>
  thread_local
  typename ref_counter<_Tp, _Alloc>::pinwall_base
  ref_counter<_Tp, _Alloc>::last(nullptr);

  template <class _Tp, template <class> class _Alloc>
  thread_local
  typename ref_counter<_Tp, _Alloc>::pinwall_base
  ref_counter<_Tp, _Alloc>::pinw(nullptr);

  template <class _Tp, template <class> class _Alloc>
  thread_local
  std::vector<typename ref_counter<_Tp, _Alloc>::counted_type*, typename ref_counter<_Tp, _Alloc>::_PtrAlloc>
  ref_counter<_Tp, _Alloc>::rmvd;

  template <class _Tp, template <class> class _Alloc>
  std::atomic<ref_data*>
  ref_counter<_Tp, _Alloc>::allrefs(nullptr);

  template <class _Tp, template <class> class _Alloc>
  thread_local
  ref_data*
  ref_counter<_Tp, _Alloc>::refdata(nullptr);

  //
  // EPOCH management scheme

  static const size_t MAXSURVIVALS = 2;

  template <class _T, class _Alloc>
  struct release_entry
  {
    typedef typename _Alloc::template rebind<size_t>::other size_t_alloc;
    typedef typename _Alloc::template rebind<_T*>::other    tp_alloc;
    typedef std::vector<size_t, size_t_alloc>               epoch_vector;
    typedef std::vector<_T*, tp_alloc>                      epoch_ptr_vector;

    epoch_vector     epochtime;
    epoch_ptr_vector epochptrs;
    size_t           survivals;

    release_entry()
    : epochtime(), epochptrs(), survivals(0)
    {}

    void swap(release_entry<_T, _Alloc>& other)
    {
      epochtime.swap(other.epochtime);
      epochptrs.swap(other.epochptrs);
    }
  };

  template <class _Tp, class _Alloc>
  struct alignas(CACHELINESZ) epoch_data
  {
    static const size_t RECOLLECTION_PERIOD = 1000; ///< collect every X operations
                                                    // \todo adjust collection period according to
                                                    //       number of threads


    typedef std::forward_list<release_entry<_Tp, _Alloc> > removal_collection;

    epoch_data<_Tp, _Alloc> * next;  ///< next thread's data
    std::atomic<size_t>       epoch; ///< odd numbers indicate busy epochs

    size_t                    collepoch;
    size_t                    stats;
    size_t                    statsctr;
    size_t                    collectctr;
    removal_collection        rmvd;

    epoch_data()
    : next(nullptr), epoch(1), collepoch(RECOLLECTION_PERIOD), stats(0), statsctr(0),
      collectctr(0), rmvd(1)
    {
      assert(!rmvd.empty());
    }
  };

  static inline
  bool active(size_t epoch)
  {
    return epoch & 1;
  }

  static inline
  bool epoch_passed(size_t freeepoch, size_t currepoch)
  {
    return (  (!active(freeepoch))    // not active at the time it was freed
           || ( currepoch > freeepoch) // or thread has advanced to new epoch
           );
  }

  template <class _Tp, class _Alloc>
  struct passed_epoch_finder
  {
    typedef epoch_data<_Tp, _Alloc>                   epoch_data_t;
    typedef release_entry<_Tp, _Alloc>                release_entry_t;
    typedef typename epoch_data_t::removal_collection removal_collection;
    typedef typename release_entry_t::epoch_vector    epoch_vector;

    typename epoch_vector::reverse_iterator curraa;

    explicit
    passed_epoch_finder(epoch_vector& vec)
    : curraa(vec.rbegin())
    {}

    bool operator()(typename removal_collection::value_type& entry)
    {
      typename epoch_vector::reverse_iterator entryzz = entry.epochtime.rend();

      return entryzz != std::mismatch(entry.epochtime.rbegin(), entryzz, curraa, epoch_passed).first;
    }
  };

  // cancels another threads transaction
  //

  template <class _Tp, class _Alloc>
  static
  void cancel_transaction(epoch_data<_Tp, _Alloc>* epochptr, size_t currepoch)
  {
#if MIXED_HTM_ATOMIC
    epochptr->epoch.compare_exchange_strong(currepoch, currepoch+1);
#else
    bool termination_success = false;

    while (!termination_success)
    {
      if (tx::begin())
      {
        if (epochptr->epoch.load(std::memory_order_relaxed) == currepoch)
        {
          epochptr->epoch.store(currepoch+1, std::memory_order_relaxed);
        }
        tx::end();
        termination_success = true;
      }
    }
#endif /* MIXED_HTM_ATOMIC */
  }

  // takes a manager to collect all epoch pointers,
  //   a range of epochs indicating when the object has been released
  //   and a pointer to a current epoch vector (note |entry| <= | curr |).
  template <class _Manager, class _Rev_iter>
  static
  void abort_delayed(_Manager& mgr, _Rev_iter entryraa, _Rev_iter entryrzz, _Rev_iter curraa)
  {
    auto vec = mgr.collect_epochptrs();
    auto raa = vec.rbegin();

    assert(std::distance(entryraa, entryrzz) <= std::distance(raa, vec.rend()));

    // we traverse the epoch vector and find all threads that have not yet
    //   made progress. if we find another thread which is active and whose
    //   epoch has not yet advanced, it will be canceled.
    while (entryraa != entryrzz)
    {
      assert(raa != vec.rend());

      size_t x = *entryraa;
      size_t y = *curraa;

      if (!epoch_passed(x, y))
      {
        cancel_transaction(*raa, *curraa);
      }

      ++entryraa; ++curraa; ++raa;
    }
  }

  template <class _Tp, template <class> class _Alloc>
  struct epoch_manager;

  template <class _Tp, template <class> class _Alloc, class _Valloc>
  struct cancel_tx_finder
  {
    typedef epoch_data<_Tp, _Alloc<_Tp> >             epoch_data_t;
    typedef typename epoch_data_t::removal_collection removal_collection;
    typedef std::vector<size_t, _Valloc>              epoch_vector;

    epoch_manager<_Tp, _Alloc>&             epochmgr;
    typename epoch_vector::reverse_iterator curraa;
#if !NDEBUG
    typename epoch_vector::reverse_iterator currzz;
#endif /*!NDEBUG*/

    explicit
    cancel_tx_finder(epoch_manager<_Tp, _Alloc>& mgr, epoch_vector& vec)
    : epochmgr(mgr), curraa(vec.rbegin())
#if !NDEBUG
      , currzz(vec.rend())
#endif /*!NDEBUG*/
    {}

    bool operator()(typename removal_collection::value_type& entry)
    {
      typedef typename epoch_vector::reverse_iterator iter_t;

      // we use a reverse iterator b/c new thread's epoch records
      //   are added at the beginning. Thus we match oldest entry with
      //   oldest entry. If an old timestamp does not have a tag
      //   for a thread, the thread must have not seen that entry.
      const iter_t              entryzz = entry.epochtime.rend();
      iter_t                    entryaa = entry.epochtime.rbegin();

      assert(std::distance(entryaa, entryzz) <= std::distance(curraa, currzz));
      std::pair<iter_t, iter_t> res = std::mismatch(entryaa, entryzz, curraa, epoch_passed);
      const bool                deletable = (entryzz == res.first);

      if (!deletable)
      {
        ++entry.survivals;

        if (entry.survivals > MAXSURVIVALS)
        {
          abort_delayed(epochmgr, entryaa, entryzz, curraa);
        }
      }

      return deletable;
    }
  };

  template <class _Tp, template <class> class _Alloc, class _Valloc>
  static
  cancel_tx_finder<_Tp, _Alloc, _Valloc>
  cancel_txepoch(epoch_manager<_Tp, _Alloc>& m, std::vector<size_t, _Valloc>& v)
  {
    return cancel_tx_finder<_Tp, _Alloc, _Valloc>(m, v);
  }


  template <class _FwdIter, class _UnaryPred>
  _FwdIter fwd_find(_FwdIter before_aa, _FwdIter zz, _UnaryPred pred)
  {
    if ( before_aa == zz ) return zz;

    _FwdIter pos = before_aa;

    ++pos;
    while ( pos != zz && !pred(*pos))
    {
      ++pos; ++before_aa;
    }

    return before_aa;
  }


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

  template <class _Alloc>
  epoch_deallocator<_Alloc> epochDeallocator(_Alloc alloc)
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
      typedef operation                 manager_kind;
      typedef alloc_base<_Tp, _Alloc>   base;
      typedef typename base::value_type value_type;

      /// no pinwall
      typedef guard<epoch_manager<value_type, _Alloc> > pinguard;

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

      void validate_transaction()
      {
        if (!htm::active(epochdata->epoch.load(std::memory_order_relaxed)))
          tx::abort<17>();
      }


      /// pins a range of pointers
      template <class _TrivialIterator>
      void txpin(_TrivialIterator, _TrivialIterator, void*)
      {
        validate_transaction();
      }

      /// unpin not needed
      void unpin(value_type*, int) {}

      //~ void node_cleanup(void (*) (value_type&), value_type&)
      //~ {
        // \todo ...
        // is the scenario in pub_scan_manager feasible under an epoch manager?
        //   T0 (epoch):   removes N0 (predecessor of N1)
        //   T1 (epoch:    removes N1 and cannot free N1
        //                ( T2 was active when it acquired N1 and will be active in the
        //                  same epoch until it finishes its operation ).
        //   T2 (active): holds a ref to N0 and pins N1 .. OK
        /* cleanupfun(n); not needed */
      //~ }

      release_entry<value_type, _Alloc<_Tp> >
      collect_epochs()
      {
        typedef scan_iterator<epoch_data<value_type, _Alloc<_Tp> > >  epoch_iterator;
        typedef epoch_data<value_type, _Alloc<_Tp> >                  epoch_data_t;
        typedef release_entry<value_type, _Alloc<_Tp> >               release_entry_t;

        epoch_data_t*    curr = allepochData.load(std::memory_order_relaxed);
        epoch_iterator   end(nullptr);
        release_entry_t  res;

        // load all epochs
        for (epoch_iterator pos(curr); pos != end; ++pos)
        {
          // \mo we load relaxed; the global barrier must have been seen before
          res.epochtime.push_back((*pos).epoch.load(std::memory_order_relaxed));
        }

        return res;
      }

      std::vector<epoch_data<value_type, _Alloc<_Tp> >*, _Alloc<epoch_data<value_type, _Alloc<_Tp> >* > >
      collect_epochptrs()
      {
        typedef epoch_data<value_type, _Alloc<_Tp> >     epoch_data_t;
        // typedef std::vector<epoch_data_t*, _Alloc<_Tp> > result_t;
        typedef std::vector<epoch_data<value_type, _Alloc<_Tp> >*, _Alloc<epoch_data<value_type, _Alloc<_Tp> >* > > result_t;
        typedef scan_iterator<epoch_data_t>              epoch_iterator;

        epoch_data_t*    curr = allepochData.load(std::memory_order_relaxed);
        epoch_iterator   end(nullptr);
        result_t         res;

        // copy all epoch pointers to the result
        for (epoch_iterator pos(curr); pos != end; ++pos)
        {
          epoch_data_t& x = *pos;

          res.push_back(&x);
        }

        return res;
      }

      void push_curr_list(typename epoch_data<value_type, _Alloc<_Tp> >::removal_collection::value_type&& epochdesc)
      {
        epochdata->rmvd.front().swap(epochdesc);
      }

      void free_unreferenced_memory()
      {
        typedef epoch_data<value_type, _Alloc<_Tp> >          epoch_data_t;
        typedef typename epoch_data_t::removal_collection     removal_collection;
        typedef typename removal_collection::iterator         removal_iterator;
        typedef release_entry<value_type, _Alloc<_Tp> >       release_entry_t;
        typedef typename release_entry_t::epoch_vector        epoch_vector;
        typedef passed_epoch_finder<value_type, _Alloc<_Tp> > passed_epoch_finder_t;

        removal_collection&    rmvd   = epochdata->rmvd;
        epoch_vector&          epoch  = rmvd.front().epochtime;
        removal_iterator       preaa  = rmvd.before_begin();
        removal_iterator const zz     = rmvd.end();

        // could use binary search (if we had a vector), but it is expected
        //   that we do not keep too many removal records around.

        // mark all unreclaimable memory as having been traversed and find
        //  position of entries where the threads should be aborted.
        removal_iterator const preabort = fwd_find(preaa, zz, cancel_txepoch(*this, epoch));

        // find the entry after which all epochs have passed and send the rest for
        //   deletion.
        removal_iterator const prebegin = fwd_find(preabort, zz, passed_epoch_finder_t(epoch));
        removal_iterator       pos      = prebegin;

        std::for_each(++pos, zz, epochDeallocator(base::get_allocator()));
        rmvd.erase_after(prebegin, zz);
      }

      /// \brief  compares epochs of releasable elements with current epochs of threads.
      ///         Frees all memory that can no longer be referenced by other threads.
      void release_memory()
      {
        assert(epochdata);

        // syncs with thread fence in the epoch counter
        std::atomic_thread_fence(std::memory_order_seq_cst);

        // collect all epochs and add current epoch to list
        push_curr_list(collect_epochs());

        // use the epochs to clear old memory
        free_unreferenced_memory();

        // add a new entry for the next epoch
        epochdata->rmvd.push_front( release_entry<value_type, _Alloc<_Tp> >() );
        ++epochdata->collectctr;
      }

      /// \brief  returns whether releasable memory blocks were held back
      ///         due to other threads being still able to reference them.
      bool has_unreleased_memory()
      {
        typedef epoch_data<value_type, _Alloc<_Tp> >      epoch_data_t;
        typedef typename epoch_data_t::removal_collection removal_collection;

        typename removal_collection::iterator aa = epochdata->rmvd.begin();
        typename removal_collection::iterator zz = epochdata->rmvd.end();

        return (  (aa != zz && (++aa) != zz)
               || (!epochdata->rmvd.front().epochptrs.empty())
               );
      }

      /// \brief  returns the number of memory blocks that are held back due to other threads
      ///         being still able to reference them.
      size_t count_unreleased_memory()
      {
         return std::accumulate( epochdata->rmvd.begin(),
                                 epochdata->rmvd.end(),
                                 0,
                                 [](size_t sz, release_entry<value_type, _Alloc<_Tp> >& entry)
                                 {
                                   return sz + entry.epochptrs.size();
                                 }
                               );
      }

      /// unpins all pointers (a noop for epoch manager)
      void unpin_all() {}

      // bool active() { return active(epochdata->epoch.load(std::memory_order_relaxed)); }

      bool active()
      {
        assert(epochdata);

        // \mo
        return htm::active(epochdata->epoch.load(std::memory_order_seq_cst));
      }


      /// increments the epoch counter
      void endop()
      {
        size_t epochval = epochdata->epoch.load(std::memory_order_relaxed);

        // exiting an operation
        if (htm::active(epochval))
        {
          epochdata->epoch.compare_exchange_strong(epochval, epochval+1, std::memory_order_relaxed, std::memory_order_relaxed);
        }
        else
        {
          // operation has been remotely terminated
        }

        // we collect every RECOLLECTION_PERIOD/2 operation cycles
        if (epochval >= epochdata->collepoch)
        {
          release_memory();

          epochdata->collepoch += epoch_data<value_type, _Alloc<_Tp> >::RECOLLECTION_PERIOD;
        }
      }

      /// \brief starts an operation and creates a data entry for the epoch managers if needed
      /// \param dummy parameter to provide consistent interface with other managers
      void beginop(size_t)
      {
        typedef epoch_data<value_type, _Alloc<_Tp> > epoch_data_t;

        if (!epochdata)
        {
          // \todo \new
          epochdata = new epoch_data_t();

          epoch_data_t* curr = allepochData.load(std::memory_order_acquire);

          do
          {
            epochdata->next = curr;
          } while (!allepochData.compare_exchange_weak(curr, epochdata, std::memory_order_release, std::memory_order_acquire)); // \mo
          return ;
        }

        // \mo this thread updates its own epoch, thus we can use relaxed
        size_t epochval = epochdata->epoch.load(std::memory_order_relaxed) + 1;

        assert(htm::active(epochval)); // entering an operation

        // \mo intra thread dependency
        epochdata->epoch.store(epochval, std::memory_order_relaxed);

        std::atomic_thread_fence(std::memory_order_seq_cst);
      }

      void deallocate(value_type* entry, int)
      {
        assert(epochdata);

        epochdata->rmvd.front().epochptrs.push_back(entry);
      }

      void stats(size_t num)
      {
        epochdata->stats += num;
        ++epochdata->statsctr;
      }

      scan_iterator<epoch_data<value_type, _Alloc<_Tp> > >
      statsbegin() { return scan_iterator<epoch_data<value_type, _Alloc<_Tp> > >(allepochData); }

      scan_iterator<epoch_data<value_type, _Alloc<_Tp> > >
      statslimit() { return scan_iterator<epoch_data<value_type, _Alloc<_Tp> > >(nullptr); }


    private:
      /// creates new storage for a Thread, unless the
      typedef std::atomic<epoch_data<value_type, _Alloc<_Tp> >*> epoch_wall_type;

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

  //
  // auxiliary functions for epochs vs non-epochs

  /// check if tx was canceled externally
  template <class _Tp, template <class> class _Alloc>
  static inline
  bool _active(epoch_manager<_Tp, _Alloc>* mgr)
  {
    // if the tx was canceled by another thread
    return (mgr->active());
  }

  /// only transactions that can be remotely canceled can be inactive
  static inline
  bool _active(void*)
  {
    return true;
  }

  //
  // StackTrack

  template <class _Tp, class _Alloc>
  struct alignas(CACHELINESZ) strack_data
  {
    typedef std::allocator_traits<_Alloc>                           _OrigAlloc_traits;
    typedef typename _OrigAlloc_traits::template rebind_alloc<_Tp*> _PpAlloc;
    typedef std::vector<_Tp*, _PpAlloc>                             removal_collection;

    std::atomic<_Tp**>        base;            ///< pinwall base
    std::atomic<_Tp**>        last;            ///< pinwall range limit
    std::atomic<size_t>       ctr;             ///< pinwall upd counter
    const size_t              pinwall_size;    ///< an estimation of the pinwall size
    strack_data*              next;            ///< next stacktrack
    size_t                    delocctr;        ///< deallocation counter
    size_t                    collection_time; ///< next dealloc timestamp
    size_t                    stats;
    size_t                    statsctr;
    size_t                    collectctr;
    removal_collection        rmvd;            ///< stores nodes to be free

    explicit
    strack_data(size_t pwsz)
    : base(nullptr), last(nullptr), ctr(0), pinwall_size(pwsz), next(nullptr),
      delocctr(0), collection_time(pwsz), stats(0), statsctr(0), collectctr(0), rmvd()
    {}
  };

  template <class _Tp, class _Alloc>
  struct strack_scanner
  {
    typedef strack_data<_Tp, _Alloc>                   strack_data_t;
    typedef typename strack_data_t::removal_collection removal_collection;

    removal_collection collect;

    void operator()(const strack_data_t& data)
    {
      // read time stamp
      // size_t lastctr = data.ctr.load(std::memory_order_relaxed);
      bool   succ = false;

      do // tx loop
      {
        // determine expected size
        _Tp**  aa = data.base.load(std::memory_order_relaxed);
        _Tp**  zz = data.last.load(std::memory_order_relaxed);

        // resize container to make sure there is enough space
        assert(aa <= zz);
        size_t expctcnt = zz-aa;

        if (expctcnt == 0) return;

        // in case we got an inconsistent read from data.base/data.load
        // \todo fix magic constant (should be szPinWall)
        if (expctcnt > 64) continue;

        // resize before to avoid the need for resizing inside a tx
        collect.reserve(collect.size() + expctcnt);

        // start transaction
        if (tx::begin())
        {
          // if time-stamp/ctr has changed -> success
          // \todo

          // if size has changed abort and restart
          if (  aa != data.base.load(std::memory_order_relaxed)
             || zz != data.last.load(std::memory_order_relaxed)
             )
          {
            tx::abort<33>();
          }

          // scan from begin to end and accumulate pointers
          // \todo we assume that a scan fits within the transaction memory
          while (aa != zz)
          {
            collect.push_back(*aa);
            ++aa;
          }

          // end transaction
          tx::end();
          succ = true;
        }
      } while (!succ);
    }

    removal_collection result() { return std::move(collect); }
  };


  template <class _Tp, template <class> class _Alloc = std::allocator >
  struct stacktrack_manager : alloc_base<_Tp, _Alloc>
  {
      typedef finegrain                 manager_kind;
      typedef alloc_base<_Tp, _Alloc>   base;
      typedef typename base::value_type value_type;
      typedef guard<stacktrack_manager<value_type, _Alloc> > pinguard;

      /// rebind to self
      template <class U>
      struct rebind { typedef stacktrack_manager<U, _Alloc> other; };

      using base::construct;
      using base::max_size;
      using base::address;

      /// constructs a publish and scan manager
      explicit
      stacktrack_manager(_Alloc<_Tp> alloc)
      : base(alloc)
      {}

      template <class _Up>
      explicit
      stacktrack_manager(stacktrack_manager<_Up, _Alloc> other)
      : base(other)
      {}

      /// constructs a publish and scan manager
      stacktrack_manager()
      : base()
      {}

    private:
      /// collects all currently published pointers from available threads
      typename strack_data<value_type, _Alloc<_Tp> >::removal_collection
      collectAllPointersFromThreads()
      {
        typedef strack_data<value_type, _Alloc<_Tp> >    strack_data;
        typedef strack_scanner<value_type, _Alloc<_Tp> > strack_scanner;
        typedef scan_iterator<strack_data>               strack_iterator;

        strack_data*    curr = allstracks.load(std::memory_order_relaxed);
        strack_iterator pos(curr);
        strack_iterator end(nullptr);
        return std::for_each(pos, end, strack_scanner()).result();
      }

    public:
      /// \brief  runs the scan and collect algorithm and frees memory that can no
      ///         longer be referenced by other threads.
      void release_memory()
      {
        typedef strack_data<value_type, _Alloc<_Tp> >    strack_data;
        typedef typename strack_data::removal_collection removal_collection;

        assert(strack);

        // syncs with thread fence in pin
        // \mo not necessary, since we sync at previous tx end or next tx start
        //~ std::atomic_thread_fence(std::memory_order_seq_cst);

        // collect and sort all currently pinned pointers
        removal_collection                    pinnedPtrs = collectAllPointersFromThreads();

        std::sort(pinnedPtrs.begin(), pinnedPtrs.end());

        // sort previously freed pointers
        removal_collection&                   rmvd = strack->rmvd;
        typename removal_collection::iterator aa   = rmvd.begin();
        typename removal_collection::iterator zz   = rmvd.end();

        // make the non pinned pointers available to be freed
        typename removal_collection::iterator pos = std::partition(aa, zz, pinned(pinnedPtrs));

        // free memory that is no longer used
        std::for_each(pos, zz, deallocator(base::get_allocator()));

        // eliminate already freed from set
        rmvd.erase(pos, zz);

        // reset collection count
        strack->delocctr = 0;
        strack->collection_time = strack->pinwall_size + threshold(pinnedPtrs.size());
        ++strack->collectctr;
      }

      /// returns the number of unreleased memory blocks
      size_t has_unreleased_memory()
      {
        assert(strack);

        return strack->rmvd.size();
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
        assert(strack);
        assert(num == 1), unused(num); // currently the allocator is limited to a single object

        strack->rmvd.push_back(obj);
        ++strack->delocctr;

        if (FREE_ALWAYS || needs_collect(strack)) release_memory();
      }

      /// \brief starts an operation and creates a data entry for the epoch managers if needed
      /// \param dummy parameter to provide consistent interface with other managers
      void beginop(size_t sz)
      {
        typedef strack_data<value_type, _Alloc<_Tp> > strack_data_t;

        if (!strack)
        {
          // \todo \new
          strack = new strack_data_t(sz);

          strack_data_t* curr = allstracks.load(std::memory_order_acquire);

          do
          {
            strack->next = curr;
          } while (!allstracks.compare_exchange_weak(curr, strack, std::memory_order_release, std::memory_order_acquire)); // \mo
        }
      }

      template <class _ForwardIterator>
      void txpin(_ForwardIterator aa, _ForwardIterator zz, void*)
      {
        htm_assert(strack);
        strack->base.store(&*aa, std::memory_order_relaxed);
        strack->last.store(&*zz, std::memory_order_relaxed);
        // ++strack->ctr;
      }

      /// unpins all entries
      void unpin_all()
      {
        htm_assert(strack);
        strack->last.store(strack->base.load(std::memory_order_relaxed), std::memory_order_relaxed);
      }

      /// ends an operation
      void endop()
      {
        unpin_all();
      }

      //
      void stats(size_t num)
      {
        strack->stats+=num;
        ++strack->statsctr;
      }

      scan_iterator<strack_data<value_type, _Alloc<_Tp> > >
      statsbegin() { return scan_iterator<strack_data<value_type, _Alloc<_Tp> > >(allstracks); }

      scan_iterator<strack_data<value_type, _Alloc<_Tp> > >
      statslimit() { return scan_iterator<strack_data<value_type, _Alloc<_Tp> > >(nullptr); }


    private:
      typedef std::atomic<strack_data<value_type, _Alloc<_Tp> >*> strack_wall_type;

      static strack_wall_type                                    allstracks;
      static thread_local strack_data<value_type, _Alloc<_Tp>>*  strack;
  };

  template <class _Tp, template <class> class _Alloc>
  typename stacktrack_manager<_Tp, _Alloc>::strack_wall_type
  stacktrack_manager<_Tp, _Alloc>::allstracks(nullptr);

  template <class _Tp, template <class> class _Alloc>
  thread_local
  strack_data<typename stacktrack_manager<_Tp, _Alloc>::value_type, _Alloc<_Tp> >*
  stacktrack_manager<_Tp, _Alloc>::strack(nullptr);

  namespace analysis
  {
    template <class Iter>
    size_t stats(Iter aa, Iter zz)
    {
      size_t sum = 0;

      for (; aa != zz; ++aa)
        sum += (*aa).stats;

      return sum;
    }

    template <class Iter>
    size_t statsctr(Iter aa, Iter zz)
    {
      size_t sum = 0;

      for (; aa != zz; ++aa)
        sum += (*aa).statsctr;

      return sum;
    }

    template <class Iter>
    size_t collectctr(Iter aa, Iter zz)
    {
      size_t sum = 0;

      for (; aa != zz; ++aa)
        sum += (*aa).collectctr;

      return sum;
    }
  }
}

#endif /* _HTMMEMORY_HPP */
