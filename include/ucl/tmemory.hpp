/// \brief     Transaction and memory managers
///
/// \details - pub_scan_manager (PS) is a memory manager that uses
///            publish and scan techniques, similar to M. Michael's hazard pointers
///            and Herlihy et al.'s repeated offender implementation.
///            This technique implies sequentially consistent reads of pointers/elements.
/// \author  Peter Pirkelbauer

#ifndef _TMEMORY_HPP

#define _TMEMORY_HPP 1

#define TXBEGIN
#define TXEND


#include <cassert>
#include <atomic>
#include <algorithm>
#include <iterator>
#include <vector>
#include <numeric>
#include <forward_list>
#include <sstream>

#include "atomicutil.hpp"
#include "bitutil.hpp"

namespace lockfree
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

      bool has_unreleased_memory() { return false; }
      size_t count_unreleased_memory() { return 0; }

      // void node_cleanup(void (*) (_Tp&), _Tp&) {}
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
  };


  //
  // Auxiliary lock-free stack implementation

  template <class _Tp, class _Alloc>
  struct alignas(CACHELINESZ) pub_scan_data
  {
      typedef std::allocator_traits<_Alloc>                          _OrigAlloc_traits;
      typedef typename _OrigAlloc_traits::template rebind_alloc<_Tp> _TpAlloc;

      typedef std::atomic<_Tp*>           pinwall_entry;
      typedef std::vector<_Tp*, _TpAlloc> removal_collection;

      pub_scan_data<_Tp, _Alloc> const * next; ///< next thread's data
      pinwall_entry* const               base; ///< base address of hazard pointers
                                               //  \todo consider to use unique pointer
      std::atomic<pinwall_entry*>        last; ///< memory managed by base
      // size_t                             maxCnt;
      size_t                             cnt;  ///< number of removed ptr since last collection
      size_t                             collection_time; ///< when the next collection should occur
      const size_t                       pinwall_size;
      removal_collection                 rmvd; ///< removed but not yet freed pointers

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

        // we confirm later anyhow, so if we do net the most recent value right
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
        assert(ofs <= 0);

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
        assert(std::find(base, zz, elem) != zz);
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

      /// returns true if this thread has unreleased memory
      bool has_unreleased_memory()
      {
        return (count_unreleased_memory() > 0);
      }

      /// returns the number of unreleased memory blocks
      size_t count_unreleased_memory()
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
        assert(num == 1); // currently the allocator is limited to a single object

        pinWall->rmvd.push_back(obj);
        ++pinWall->cnt;

        if (FREE_ALWAYS || pinWall->needs_collect()) release_memory();
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


        // set txsize
        txstarted = false;
        txsize = TXMAGICNUMBER;
        while (txsize > 1 || txstarted)
        {
          if (TXBEGIN)
          {
            // regular transaction
            txstarted = true;
          }
          else
          {
            // conflict handling + retry
            txsize = txsize / 2;
          }
        }
      }

      /// \brief pins an atomic pointer
      template <std::memory_order mo = std::memory_order_seq_cst>
      value_type*
      pin(std::atomic<value_type*>& elem)
      {
        if (txsize == 1)
        {
          // do the real stuff
          assert(pinWall);
          return pinWall->template pin_aux<mo>(elem);
        }
        else
        {
          // to the tx stuff
          //   i.e., simply read
        }
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

      /// node cleanup invokes a specified cleanup function
      ///   to clean up a node that can be freed.
      ///   In non GCed environments links to next nodes need to be set to null,
      ///   otherwise some delayed thread may access the next element, which
      ///   may have been freed previously, since the next element is in noone's
      ///   published list.
      ///   T0: removes N0 (predecessor of N1)
      ///   T1: removes N1 and frees N1 (it is on no published list)
      ///   T2: holds a ref to N0 and pins N1 ... oops
      //~ void node_cleanup(void (cleanupfun) (value_type&), value_type& n)
      //~ {
        //~ cleanupfun(n);
      //~ }

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
}

#endif /* _PMEMORY_HPP */
