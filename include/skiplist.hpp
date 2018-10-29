/// \file   skiplist.hpp
/// \author Nick Dzugan, Peter Pirkelbauer
///
/// \brief  Lock-free and finegrain-locking skiplist implementations.
/// \details
///         - namespace locking: fine-grained locking skiplist
///         - namespace lockfree: lock-free skiplist
///
///         The implementations are based on Herlihy and Shavit: The Art of
///         Multiprocessor Programming, 2012.

#ifndef _SKIPLIST_HPP
#define _SKIPLIST_HPP 1

#include <memory>
#include <iterator>
#include <mutex>
#include <algorithm>
#include <array>
#include <random>
#include <climits>

#include "pmemory.hpp"
#include "bitutil.hpp"

/// \private
namespace aux
{
  /// \private
  thread_local std::ranlux48 rnd(time(0));

  /// \private
  static inline
  int bsr32(uint32_t s)
  {
    if (s == 0) return 0;

    static const int nobits = (sizeof(unsigned int) * CHAR_BIT);

    // note: use __builtin_clzl, and __builtin_clzll for long and long long
    return nobits - __builtin_clz(s);
  }

  /// \private
  constexpr
  size_t _next_power2(size_t N, size_t X)
  {
    return (X > 16) ? N : _next_power2(N | (N >> X), X << 1);
  }

  /// \private
  constexpr
  size_t next_power2(size_t N)
  {
    return _next_power2(N-1, 1) + 1;
  }

/*
  template <size_t M>
  static inline
  size_t log_rand()
  {
    const size_t X = next_power2(M);
    const size_t mask = (uint64_t(1) << (X)) - 1;
    const size_t r = X - bsr32(rnd() & mask);

    return (r < M) ? r : log_rand<M>();
  }
*/

  /// \private
  template <size_t M>
  static inline
  size_t log_rand()
  {
    const size_t X = next_power2(M);
    const size_t mask = (uint64_t(1) << (X)) - 1;
    const size_t r = (X - bsr32(rnd() & mask)) / 3;

    return (r < M) ? r : log_rand<M>();
  }

  /// \private
  template <class _NodeTp>
  class skiplist_iterator : std::iterator< std::forward_iterator_tag, typename _NodeTp::value_type>
  {
    typedef std::iterator< std::forward_iterator_tag, typename _NodeTp::value_type> base;

    public:
      typedef typename base::value_type        value_type;
      typedef typename base::reference         reference;
      typedef typename base::pointer           pointer;
      typedef typename base::difference_type   difference_type;
      typedef typename base::iterator_category iterator_category;

      explicit
      skiplist_iterator(_NodeTp* node)
      : skipNode(node)
      {}

      /// \private
      template <class U>
      friend
      bool operator==(const skiplist_iterator<U>& lhs, const skiplist_iterator<U>& rhs);

      /// \private
      template <class U>
      friend
      bool operator!=(const skiplist_iterator<U>& lhs, const skiplist_iterator<U>& rhs);

      reference operator*() const;
      skiplist_iterator<_NodeTp>& operator++();
      skiplist_iterator<_NodeTp>  operator++(int);

    private:
      _NodeTp* skipNode;
  };

  /// \private
  template <class T>
  bool operator==(const skiplist_iterator<T>& lhs, const skiplist_iterator<T>& rhs)
  {
    return (lhs.skipNode == rhs.skipNode);
  }

  /// \private
  template <class T>
  bool operator!=(const skiplist_iterator<T>& lhs, const skiplist_iterator<T>& rhs)
  {
    return (lhs.skipNode != rhs.skipNode);
  }

  template <class T>
  typename skiplist_iterator<T>::reference
  skiplist_iterator<T>::operator*() const
  {
    return skipNode->nodeValue;
  }

  template <class T>
  skiplist_iterator<T>& skiplist_iterator<T> :: operator++()
  {
    skipNode = skipNode->next[0].load(std::memory_order_relaxed); // \mo
    return *this;
  }

  template <class T>
  skiplist_iterator<T> skiplist_iterator<T> :: operator++(int)
  {
    skiplist_iterator<T> original(*this);
    skipNode = skipNode->next[0].load(std::memory_order_relaxed); // \mo
    return original;
  }
}

namespace locking
{
  /// \private
  template <class T, class _Alloc>
  class skiplist_node
  {
    typedef std::atomic<skiplist_node<T, _Alloc>*> node_ptr;

    typedef std::allocator_traits<_Alloc>                                _OrigAlloc_traits;
    typedef typename _OrigAlloc_traits::template rebind_alloc<node_ptr>  _NodePtr_alloc_type;

    public:
      typedef T value_type;

      // std::atomic<bool>          isdead;
      std::recursive_mutex       myLock;      ///< per element lock
                                              ///  since a thread can acquire the same
                                              ///  at multiple levels of the skiplist,
                                              ///  we use a recursive lock.
                                              //  \todo since locks are acquired through
                                              //        direct nesting, we could optimize
                                              //        and use a simpler mutex
      T                          nodeValue;   ///< user data
      std::atomic<bool>          marked;      ///< marked for deletion
      std::atomic<bool>          fullyLinked; ///< node has been set up and is valid
      node_ptr*                  next;        ///< linkage levels; the number of
                                              ///  linkage levels is computed randomly
      const size_t               noLevels;    ///< number of linkage levels

      // remove auto-generated ctors
      skiplist_node() = delete;
      skiplist_node(const skiplist_node&) = delete;
      skiplist_node& operator=(const skiplist_node&) = delete;

      /// Constructs a new skiplist node for an entry newentry; the node will
      /// participate on levelsLinked levels
      skiplist_node(const T& newentry, int levelsLinked)
      : // isdead(false),
        myLock(), nodeValue(newentry),
        marked(false), fullyLinked(false),
        next(nullptr), noLevels(levelsLinked)
      {
        assert(levelsLinked > 0);
        _NodePtr_alloc_type alloc;

        next = alloc.allocate(noLevels);
      }

      ~skiplist_node()
      {
        _NodePtr_alloc_type alloc;

        alloc.deallocate(next, noLevels);
      }

      /// locks node for manipulation
      void lock()
      {
        // assert(!isdead.load(std::memory_order_relaxed));
        myLock.lock();
      }

      /// unlocks node
      void unlock()
      {
        myLock.unlock();
      }

      /// number of linkage levels
      size_t levels() const { return noLevels; }
  };

  /// skiplist implementation
  /// \brief   A skiplist based set implementation that uses fine-grained locking.
  /// \tparam _Tp element type
  /// \tparam _Compare less-than operation
  /// \tparam _Alloc allocator; the default allocator uses a hazard pointer scheme
  /// \tparam _MAXLVL maximum number of levels in skiplist ( 0 < _MAXLVL <= 32 )
  template < class  _Tp,
             class  _Compare = std::less<_Tp>,
             class  _Alloc = lockfree::epoch_manager<_Tp, std::allocator>,
             size_t _MAXLVL = 32
           >
  class skiplist
  {
    public:
      static const int MAXLEVEL  = _MAXLVL;
      static const int MAX1      = MAXLEVEL+1;
      static const int NUMLOCALS = 0; // \pp was 2;
      static const int SZPINWALL = 2*MAX1 + NUMLOCALS;

      typedef skiplist_node<_Tp, typename _Alloc::base_allocator> node_type;

      // \todo add boost concept checks

      typedef typename _Alloc::value_type                                   _Alloc_value_type;
      typedef std::allocator_traits<_Alloc>                                 _OrigAlloc_traits;
      typedef typename _OrigAlloc_traits::template rebind_alloc<_Tp>        _Tp_alloc_type;
      typedef typename _OrigAlloc_traits::template rebind_alloc<node_type>  _Node_alloc_type;
      typedef std::allocator_traits<_Tp_alloc_type>                         _Alloc_traits;
      typedef typename _Node_alloc_type::pinguard                           PinGuard;

    public:
      typedef _Tp                                             value_type;
      typedef _Tp                                             key_type;

      typedef _Compare                                        value_compare;
      typedef _Compare                                        key_compare;

      typedef typename _Alloc_traits::pointer                 pointer;
      typedef typename _Alloc_traits::const_pointer           const_pointer;
      typedef typename _Alloc_traits::value_type&             reference;
      typedef const typename _Alloc_traits::value_type&       const_reference;

      typedef aux::skiplist_iterator<node_type>               iterator;
      typedef iterator                                        const_iterator; // \todo make real const iterator
      typedef typename iterator::difference_type              difference_type;

      typedef size_t                                          size_type;
      typedef _Alloc                                          allocator_type;

      explicit
      skiplist( const key_compare& comparator = key_compare(),
                const allocator_type& alloc = allocator_type()
              )
      : start(value_type(), MAXLEVEL),
        limit(value_type(), MAXLEVEL),
        lessThan(comparator),
        nodeAlloc(alloc)
      {
        // must be at least one level
        static_assert(MAXLEVEL >  0,  "MAXLEVEL > 0 required");

        // we use bsr32 in log_rand, hence the max level is 32
        static_assert(MAXLEVEL <= 32, "MAXLEVEL <= 32 required");

        // fine-grained memory management techniques do not work, when we use
        //   optimistic read-through operations.
        // \todo explain in documentation
        static_assert(!std::is_same<typename _Alloc::manager_kind, lockfree::finegrain>::value, "fine grained memory manager do not work with this skiplist implementation (try lockfree::skiplist instead).");

        for (int i = 0; i < MAXLEVEL; ++i)
        {
          start.next[i].store(&limit, std::memory_order_relaxed);  // \mo
          limit.next[i].store(nullptr, std::memory_order_relaxed); // \mo
        }

        start.fullyLinked.store(true, std::memory_order_relaxed); // \mo
        limit.fullyLinked.store(true, std::memory_order_relaxed); // \mo

        // \todo add barrier to prevent reordering w/ code outside ctor?
      }

      /// \brief inserts a new entry into the list
      /// \return true, if the entry was inserted; false, otherwise.
      bool insert(const value_type& newentry);

      /// \brief removes entry from the list.
      /// \return true, if entry was in the list; false otherwise.
      bool erase(const value_type& entry);

      /// \brief tests whether entry is in the skiplist
      bool contains(const value_type& entry);

      /// \brief  returns the allocator of this skiplist
      _Node_alloc_type getAllocator() const noexcept
      {
        return nodeAlloc;
      }

      /// \brief returns an iterator to the first element
      /// \details the function is only safe as long as the skip-list is in a quiescent state
      iterator qbegin();

      /// \brief returns an iterator to the last element
      /// \details the function is only safe as long as the skip-list is in a quiescent state
      iterator qend();

#if TODO
      const_iterator qcbegin();
      const_iterator qcend();
#endif /* TODO */

    private:
      /// Auxiliary Method to find the location of a given element
      int find (const value_type& x, std::array<node_type*, MAX1>& preds, std::array<node_type*, MAX1>& succs);

      node_type           start;
      node_type           limit;
      _Compare            lessThan;
      _Node_alloc_type    nodeAlloc;

      void lock( node_type* n )
      {
        n->lock();
      }
  };


  /// \private
  template <int MAX, class _NodeType>
  static
  void unlockLevels(std::array<_NodeType*, MAX>& preds, int lb, int ub)
  {
    assert(lb >= 0 && ub >= lb && ub < MAX);

    for (int level = lb; level < ub; ++level)
    {
      preds[level]->unlock();
    }
  }

  template <class _Tp, class _Compare, class _Alloc, size_t _MAXLVL>
  bool skiplist<_Tp, _Compare, _Alloc, _MAXLVL> :: insert(const value_type& x)
  {
    const int maxLevel = aux::log_rand<MAXLEVEL>()+1;
    assert(maxLevel > 0 && maxLevel <= MAXLEVEL);

    std::array<node_type*, MAX1> preds;
    std::array<node_type*, MAX1> succs;

    while (true)
    {
      PinGuard  guard(nodeAlloc, SZPINWALL); // uses RAII to free alloc's pinwall
      int       lFound = find(x, preds, succs);

      if (lFound != -1)
      {
        // We found a node with the same key
        node_type* nodeFound = succs[lFound];

        // \mo
        // sequentially consistent, since we did not acquire a lock
        //   syncs with marked.store in erase
        if (!nodeFound->marked.load())
        {
          // wait for it to become fully linked
          // \mo relaxed b/c eventually the new value becomes available
          while (!nodeFound->fullyLinked.load(std::memory_order_relaxed)) {}

          return false;
        }

        // The node is marked for deletion
        //   -> repeat until we can insert
        continue;
      }

      int        highestLocked = -1;
      node_type* pred = nullptr;
      node_type* succ = nullptr;
      bool       valid = true;

      for (int level = 0; valid && (level < maxLevel); ++level)
      {
        pred = preds[level];
        succ = succs[level];

        // here we need a reentrant lock, as it is possible to lock the same
        //   node on multiple levels.
        // \todo try to use hierarchical info to avoid secondary lock acquisitions
        //       and replace lock with simpler mutex
        lock(pred);
        highestLocked = level;

        valid = (  !pred->marked.load(std::memory_order_relaxed) // \mo through lock acq
                && !succ->marked.load(std::memory_order_seq_cst) // \mo there is no other
                                                                 //     ordering relationship with
                                                                 //     succ, thus we rely on
                                                                 //     sequential consistency
                && pred->next[level].load(std::memory_order_relaxed) == succ // \mo through lock acq
                /* && succ->prev == pred; */
                );
      }

      if (!valid)
      {
        unlockLevels<MAX1>(preds, 0, highestLocked+1);
        continue;
      }

      // allocate & construct new node
      _Node_alloc_type alloc   = getAllocator();
      node_type*       newNode = alloc.allocate(1);

      alloc.construct(newNode, x, maxLevel);

      // interlink new node
      for (int level = 0; level < maxLevel; ++level)
      {
        // \mo a node may not be visible before it gets fully linked
        //     thus we can utilize relaxed memory operations
        //     the synchronization point is deferred to setting fullyLinked
        //     to true below.
        newNode->next[level].store(succs[level], std::memory_order_relaxed);
      }

      // link node into skiplist
      for (int level = 0; level < maxLevel; ++level)
      {
        // \mo release to guarantee that the nodes linkage to successors
        //     and the linkages on the lower levels are visible
        preds[level]->next[level].store(newNode, std::memory_order_release);
      }

      // may need to move.
      newNode->fullyLinked.store(true, std::memory_order_seq_cst); // \mo syncs with loading fully_linked in erase

      // unlock all
      unlockLevels<MAX1>(preds, 0, highestLocked+1);
      return true;
    }
  }


  template <class _Tp, class _Compare, class _Alloc, size_t _MAXLVL>
  int
  skiplist<_Tp, _Compare, _Alloc, _MAXLVL>
  ::find(const value_type& x, std::array<node_type*, MAX1>& preds, std::array<node_type*, MAX1>& succs)
  {
    int              levelFound = -1;
    _Node_alloc_type alloc      = getAllocator();

    // \note: While &start can never be freed (w/o freeing the skiplist),
    //        we want to treat it uniformly with other pointers (can be unpinned).
    node_type*       pred       = &start;

    for (int level = MAXLEVEL-1; level >= 0; --level)
    {
      alloc.pin_addr(*pred);

      node_type* curr = alloc.template pin<std::memory_order_consume>(pred->next[level]);
      // assert( curr != nullptr ); // \note does not hold under publish and collect managers

      while (curr != nullptr && curr != &limit && lessThan(curr->nodeValue, x))
      {
        alloc.unpin(pred, -1);
        pred = curr;
        curr = alloc.template pin<std::memory_order_consume>(pred->next[level]);
      }

#if ERRONEOUS_CODE
      // Under P&S manager, the algorithm may find a null-ptr (node marked
      //   to be freed), in which case we restart the algorithm.
      // \todo consider more fine-grained ways of restarting the algorithm
      //       or removing support for finegrain memory management techniques.
      if (curr == nullptr)
      {
        // manually end operation, b/c slots will be needed in recursive call
        return find(x, preds, succs);
      }
#endif /* ERRONEOUS_CODE */

      if (  levelFound == -1
         && curr != &limit
         && !lessThan(x, curr->nodeValue)
         )
      {
        levelFound = level;
      }

      preds[level] = pred;
      succs[level] = curr;
    }

    return levelFound;
  }

#if ERRONEOUS_CODE
  template <class _Tp>
  static
  void unlinkNodeCleanup(skiplist_node<_Tp>& node)
  {
    for (int lv = node.levels()-1; lv >= 0; --lv)
    {
      node.next[lv].store(nullptr, std::memory_order_relaxed);
    }
  }
#endif /* ERRONEOUS_CODE */

  template <class _Tp, class _Compare, class _Alloc, size_t _MAXLVL>
  bool skiplist<_Tp, _Compare, _Alloc, _MAXLVL> :: erase(const value_type& x)
  {
    bool                         isMarked = false;
    int                          maxLevel = -1;
    std::array<node_type*, MAX1> preds;
    std::array<node_type*, MAX1> succs;

    // \todo can we replace the calls to unlockLevels with an RAII guard?
    while (true)
    {
      PinGuard   guard(nodeAlloc, SZPINWALL); // uses RAII to free alloc's pinwall
      int        lFound = find(x, preds, succs);
      node_type* victim = nullptr;

      if (lFound != -1) //if is found
      {
        victim = succs[lFound];
      }

      // \mo At this point the lock for victim is not yet acquired.
      // Thus we rely on sequential consistency to order operations on
      // fullyLinked and marked
      if (  isMarked
         || (  victim                     // if node is found in the list
            && victim->fullyLinked.load()
            && (victim->levels()-1) == static_cast<size_t>(lFound)
            && !victim->marked.load()
            ) )
      {
        if (!isMarked)
        {
          maxLevel = victim->levels();
          assert(maxLevel > 0);

          lock(victim); // \todo may not be able to lock the victim. may have to rely on marking.

          // \mo
          // relaxed load, b/c if we load here we have acquired the lock
          //   previously. Another thread setting marked would have
          //   released the lock earlier.
          if (victim->marked.load(std::memory_order_relaxed))
          {
            victim->unlock();
            return false;
          }

          victim->marked.store(true, std::memory_order_seq_cst); // \mo syncs with load in condition above
          isMarked = true;
        }

        int highestLocked = -1;

        // \todo Block could use RAII to unlock the locked nodes instead of
        //       calling unlock-levels explicitly
        {
          node_type* pred = nullptr;
          bool       valid = true;

          for (int level = 0; valid && (level < maxLevel); ++level)
          {
            pred = preds[level];
            lock(pred);
            highestLocked = level;

            // \mo
            // relaxed after lock acq on pred
            valid = (  !pred->marked.load(std::memory_order_relaxed)
                    && pred->next[level].load(std::memory_order_relaxed) == victim
                    );
          }

          if (!valid)
          {
            unlockLevels<MAX1>(preds, 0, highestLocked+1);
            continue;
          }

          for (int level = maxLevel-1; level >= 0; --level)
          {
            // \mo
            // relaxed, since we are holding the lock
            node_type* succ = victim->next[level].load(std::memory_order_relaxed);

            // \mo
            // release to guarantee that a reading thread sees unlinking the node
            // in the proper sequence
            preds[level]->next[level].store(succ, std::memory_order_release);
          }

          unlockLevels<MAX1>(preds, 0, highestLocked+1); // \todo move up

          _Node_alloc_type alloc = getAllocator();


#if ERRONEOUS_CODE
          // Unlink victim for publish-and-scan managements
          alloc.node_cleanup(unlinkNodeCleanup<value_type>, *victim);
#endif /* ERRONEOUS_CODE */

          // free victim
          victim->unlock();

          /* alloc.destroy(victim); no longer needed */
          alloc.deallocate(victim, 1);
          return true;
        }
      }
      else // \todo can the else be removed and the return just be put at the end of the block?
      {
        return false;
      }
    }

    assert(false); // cannot be reached
  }


  template <class _Tp, class _Compare, class _Alloc, size_t _MAXLVL>
  bool skiplist<_Tp, _Compare, _Alloc, _MAXLVL> :: contains(const value_type& x)
  {
    std::array<node_type*, MAX1> preds;
    std::array<node_type*, MAX1> succs;
    PinGuard                     guard(nodeAlloc, SZPINWALL); // uses RAII to free alloc's pinwall
    const int                    lFound = find(x, preds, succs);

    // \mo
    // sequentially consistent, since we do not acquire any locks
    return (  lFound != -1
           && succs[lFound]->fullyLinked.load()
           && !succs[lFound]->marked.load()
           );
  }

#if TODO
  // \pp \todo needs to be implemented in terms of find...
  template <class _Tp, class _Compare, class _Alloc, size_t _MAXLVL>
  skiplist_iterator<_Tp> skiplist<_Tp, _Compare, _Alloc, _MAXLVL> (const value_type& x)
  {
      node_type* newIterator = &start;

      while (newIterator != &limit)
      {
          if (newIterator->nodeValue == x)
              return iterator(newIterator);
          else
              newIterator = newIterator->next[0].load(); // \mo
      }
      return end();
  }
#endif /* TODO */

  template <class _Tp, class _Compare, class _Alloc, size_t _MAXLVL>
  typename skiplist<_Tp, _Compare, _Alloc, _MAXLVL>::iterator
  skiplist<_Tp, _Compare, _Alloc, _MAXLVL> :: qbegin()
  {
    return iterator(start.next[0].load());
  }

  template <class _Tp, class _Compare, class _Alloc, size_t _MAXLVL>
  typename skiplist<_Tp, _Compare, _Alloc, _MAXLVL>::iterator
  skiplist<_Tp, _Compare, _Alloc, _MAXLVL> :: qend()
  {
    return iterator(&limit);
  }
} // namespace locking

namespace lockfree
{
  /// \private
  template <class _Tp, class _Alloc>
  class skiplist_node
  {
      typedef uab::MarkablePointer<skiplist_node<_Tp, _Alloc>>            node_ptr;
      typedef std::allocator_traits<_Alloc>                               _OrigAlloc_traits;
      typedef typename _OrigAlloc_traits::template rebind_alloc<node_ptr> _NodePtr_alloc_type;

    public:
      typedef _Tp value_type;

      // std::atomic <size_t>       isdead;        //is is dead necessary?
      _Tp          node_value;   ///< user data
      const size_t num_levels;   ///< number of linkage levels
      node_ptr*    next;         ///< pointer of array that connext to the successors (the depth is computed randomly)

      /// Constructs a new skiplist node for an entry new_entry; the node will
      /// participate on levels_linked levels
      skiplist_node(const _Tp& new_entry, int levels_linked)
      : node_value(new_entry),num_levels(levels_linked), next(nullptr)
      {
        assert(levels_linked > 0);
        _NodePtr_alloc_type alloc;

        next = alloc.allocate(num_levels);
      }

      ~skiplist_node()
      {
        _NodePtr_alloc_type alloc;

        alloc.deallocate(next, num_levels);
      }

      /// number of linkage levels
      size_t levels() const { return num_levels; }
  };

  /// \brief a lock-free skiplist. Depending on the memory manager,
  ///        the skiplist may also be nonblocking.
  /// \tparam _Tp element type
  /// \tparam _Compare less-than operation
  /// \tparam _Alloc allocator;
  /// \tparam _MAXLVL maximum number of levels in skiplist (0 < _MAXLVL <= 32)
  template < class  _Tp,
             class  _Compare = std::less<_Tp>,
             class  _Alloc = lockfree::epoch_manager<_Tp, std::allocator>,
             size_t _MAXLVL = 32
           >
  class skiplist
  {
      typedef typename uab::MarkablePointer<skiplist_node<_Tp,_Alloc>*>::mark_type mark_type;

    public:
      static const int MAXLEVEL  = _MAXLVL;

    private:
      static const int MAX1      = MAXLEVEL+1;
      static const int NUMLOCALS = 0; // \pp was 2
      static const int SZPINWALL = 3*MAX1 + NUMLOCALS;

      typedef skiplist_node<_Tp, typename _Alloc::base_allocator>           node_type;
      typedef typename _Alloc::value_type                                   _Alloc_value_type;
      typedef std::allocator_traits<_Alloc>                                 _OrigAlloc_traits;
      typedef typename _OrigAlloc_traits::template rebind_alloc<_Tp>        _Tp_alloc_type;
      typedef typename _OrigAlloc_traits::template rebind_alloc<node_type>  _Node_alloc_type;
      typedef std::allocator_traits<_Tp_alloc_type>                         _Alloc_traits;
      typedef typename _Node_alloc_type::pinguard                           MemGuard;

    public:
      typedef _Tp                                                           value_type;
      typedef _Tp                                                           key_type;

      typedef _Alloc                                                        allocator_type;
      typedef aux::skiplist_iterator<node_type>                             iterator;
      typedef iterator                                                      const_iterator; // \todo make real const iterator
      typedef typename iterator::difference_type                            difference_type;

      typedef _Compare                                                      value_compare;
      typedef _Compare                                                      key_compare;

      typedef typename _Alloc_traits::pointer           pointer;
      typedef typename _Alloc_traits::const_pointer     const_pointer;
      typedef typename _Alloc_traits::value_type&       reference;
      typedef const typename _Alloc_traits::value_type& const_reference;

      explicit
      skiplist( const key_compare& comparator = key_compare(),
                const allocator_type& alloc = allocator_type()
              )
      : start(value_type(), MAXLEVEL),
        limit(value_type(), MAXLEVEL),
        lessThan(comparator),
        nodeAlloc(alloc)
      {
        static_assert(MAXLEVEL >  0,  "MAXLEVEL > 0 required");
        static_assert(MAXLEVEL <= 32, "MAXLEVEL <= 32 required");

        for (int i = 0; i < MAXLEVEL; ++i)
        {
          start.next[i].store(&limit, std::memory_order_relaxed);  // \mo
          limit.next[i].store(nullptr, std::memory_order_relaxed); // \mo
        }

        // \todo add barrier?
      }

      /// \brief inserts a new entry into the list
      /// \return true, if the entry was inserted; false, otherwise.
      bool insert(const value_type& new_entry);

      /// \brief removes entry from the list.
      /// \return true, if entry was in the list; false otherwise.
      bool erase(const value_type& entry);

      // skiplist_node<_Tp>* pop();

      /// \brief size of the skiplist
      /// \details the function is only safe as long as the skip-list is in a quiescent state
      size_t qsize() /* \todo const */;

      /// \brief tests whether entry is in the skiplist
      bool contains(const value_type& entry);

      /// \brief  returns the allocator of this skiplist
      _Node_alloc_type getAllocator() const noexcept
      {
        return nodeAlloc;
      }

      /// \brief returns an iterator to the first element
      /// \details the function is only safe as long as the skip-list is in a quiescent state
      iterator qbegin();

      /// \brief returns an iterator to the last element
      /// \details the function is only safe as long as the skip-list is in a quiescent state
      iterator qend();

#if TODO
      const_iterator qcbegin();
      const_iterator qcend();
#endif /* TODO */

    private:
      /// Auxiliary Method to find the location of a given element
      std::pair<node_type*, node_type*>
      find (const value_type& x, std::array<node_type*, MAX1>& preds, std::array<node_type*, MAX1>& succs);

      node_type        start;
      node_type        limit;
      _Compare         lessThan;
      _Node_alloc_type nodeAlloc;
  };


  template <class _Tp, class _Compare, class _Alloc, size_t _MAXLVL>
  bool skiplist<_Tp, _Compare, _Alloc, _MAXLVL> :: insert(const value_type& x)
  {
    const int maxLevel = aux::log_rand<MAXLEVEL>()+1;
    assert(maxLevel > 0 && maxLevel <= MAXLEVEL);

    uintptr_t                    _false = false;
    std::array<node_type*, MAX1> preds;
    std::array<node_type*, MAX1> succs;

    int bottom_level = 0;

    while (true)
    {
      MemGuard                          guard(nodeAlloc, SZPINWALL); // uses RAII to free alloc's pinwall
      std::pair<node_type*, node_type*> loc = find(x, preds, succs);

      if (loc.first != nullptr)   // x is already in the list
        return false;

      _Node_alloc_type alloc   = getAllocator();
      node_type*       new_node = alloc.allocate(1);

      alloc.construct(new_node, x, maxLevel);

      // use relaxed to set up node
      // \mo setup uses relaxed as the node is not yet visible
      //     uses release once the bottom level is linked in
      for (int level = bottom_level; level < maxLevel; ++level)
      {
        new_node->next[level].store(succs[level], std::memory_order_relaxed); // \mo
      }

      node_type* pred = preds[bottom_level];
      node_type* succ = succs[bottom_level];
      _false = false;

      new_node->next[bottom_level].store(succs[bottom_level], std::memory_order_relaxed);

      if (!pred->next[bottom_level].compare_exchange_strong( succ,
                                                             new_node,
                                                             _false,
                                                             false,
                                                             std::memory_order_release,
                                                             std::memory_order_relaxed
                                                           ))
      {
        // if the bottom level has changed, restart
        continue;
      }

      // \todo why not start at bottom_level + 1?
      for (int level = bottom_level; level < maxLevel; ++level)
      {
        while (true)
        {
          pred = preds[level];
          succ = succs[level];
          _false = false;

          // \1 release b/c we publish our data
          // \2 relaxed b/c the operation did not succeed and we have to reread preds and succs
          if (pred->next[level].compare_exchange_strong( succ,
                                                         new_node,
                                                         _false,
                                                         false,
                                                         std::memory_order_release,
                                                         std::memory_order_relaxed
                                                       ))
          {
            break;
          }

          // pinned memory must be cleared before we can reexecute find
          //   (needed in order not to exceed pointer slots)
          alloc.unpin_all();
          find(x, preds, succs);

          // \todo do not we have to set the successor nodes, if something has changed?
        }
      }

      return true;
    }
  }

  template <class _Tp, class _Compare, class _Alloc, size_t _MAXLVL>
  auto
  skiplist<_Tp, _Compare, _Alloc, _MAXLVL>
  ::find(const value_type& x, std::array<node_type*, MAX1>& preds, std::array<node_type*, MAX1>& succs)
  -> std::pair<node_type*, node_type*>
  {
    typedef std::pair<node_type*, mark_type> ptr_state;

    static const mark_type F = false;
    static const int       bottom_level = 0;

    _Node_alloc_type alloc = getAllocator();

    retry:
    {
      node_type* pred = &start;
      ptr_state  curr;
      // keep track of save pointers
      // \todo separate out into allocation manager's class

      for (int level = MAXLEVEL-1; level >= bottom_level; --level)
      {
        size_t C = -1;
        size_t P = -2;

        alloc.pin_addr(*pred);

        // \mo consume b/c we will access data through the loaded pointer
        curr = alloc.template pin<std::memory_order_consume>(pred->next[level]);

        while (curr.first != &limit)
        {

          // \mo consume b/c we may access data through the loaded pointer
          ptr_state succ = alloc.template pin<std::memory_order_consume>(curr.first->next[level]);

          // if this node is marked, we want to remove it
          while (succ.second)
          {
            // we expect pred unmarked, hence falseval
            mark_type falseval = F;

            // \mo \1 release to make available succ's data
            //     \2 relaxed b/c op has no effect and we are going to restart
            bool      snip = pred->next[level].compare_exchange_strong( curr.first,
                                                                        succ.first,
                                                                        falseval,
                                                                        F,
                                                                        std::memory_order_release,
                                                                        std::memory_order_relaxed
                                                                      );

            // if pred's next node has changed (either different successor node, or marked bit)
            //   we restart.
            if (!snip)
            {
              alloc.unpin_all();
              goto retry;
            }

            alloc.unpin(curr.first, C);
            curr = succ;

            // reached the end
            if (curr.first == &limit) break;

            // \mo see previous comment on loading succ
            succ = alloc.template pin<std::memory_order_consume>(curr.first->next[level]);
          }

          if ((curr.first == &limit) || !lessThan(x, curr.first->node_value)) break;

          alloc.unpin(pred, P);

          // \todo Here the p&s implementation leaks into the list.
          //       find better ways for an efficient p&s implementation!
          std::swap(P, C);

          pred = curr.first;
          curr.first = succ.first;
        }

        preds[level] = pred;
        succs[level] = curr.first;
      }

      if ((curr.first != &limit) && !lessThan(curr.first->node_value, x))
      {
        return std::make_pair(pred, curr.first);
      }

      return std::make_pair(nullptr, nullptr);
    }
  }

  template <class _Tp, class _Compare, class _Alloc, size_t _MAXLVL>
  bool skiplist<_Tp, _Compare, _Alloc, _MAXLVL> :: erase(const value_type& x)
  {
    typedef std::pair<node_type*, mark_type> ptr_state;

    static const int             bottom_level = 0;

    _Node_alloc_type             alloc = getAllocator();
    std::array<node_type*, MAX1> preds;
    std::array<node_type*, MAX1> succs;

    while (true)
    {
      MemGuard                          guard(alloc, SZPINWALL); // uses RAII to free alloc's pinwall
      std::pair<node_type*, node_type*> loc = find(x, preds, succs);

      // if not found
      if (loc.first == nullptr) return false;

      // found the element
      node_type* victim = succs[bottom_level];

      // mark from top to bottom (except for the lowest level)
      for (int level = victim->levels()-1; level > bottom_level; --level)
      {
        // loop until level marked
        for (;;)
        {
          // \mo relaxed b/c if someone else has marked it, we will see in the next loop
          ptr_state succ = victim->next[level].state(std::memory_order_relaxed);

          if (succ.second) break;

          victim->next[level].mark_strong( succ.first,
                                           true,
                                           std::memory_order_release,
                                           std::memory_order_relaxed
                                         );
        }
      }

      // \mo relaxed, b/c we will execute the loop until a change has been made
      ptr_state succ = victim->next[bottom_level].state(std::memory_order_relaxed);

      while (true)
      {
        // \mo
        // \1 release guarantees that all previous mark bits are found
        // \2 relaxed, b/c someone else modified the node
        bool i_marked_it = victim->next[bottom_level].mark_strong( succ.first,
                                                                   true,
                                                                   std::memory_order_release,
                                                                   std::memory_order_relaxed
                                                                 );

        if (i_marked_it)
        {
          // call find again to eventually unlink all nodes
          alloc.unpin_all();
          find(x, preds, succs);
          return true;
        }

        // Check state of node
        // \mo \todo is there any dependent data to use acquire?
        succ = victim->next[bottom_level].state(std::memory_order_consume);

        // This node has been deleted by someone else
        if (succ.second)
          return false;

        // Someone else has changed the successor, so we re-execute the loop
      }
    }
  }

#if 0
  template <class _Tp, class _Compare, class _Alloc, size_t _MAXLVL>
  skiplist_node<_Tp>* skiplist<_Tp, _Compare, _Alloc, _MAXLVL> :: pop()
  {
    // \pir \todo add memory management

    int bottom_level = 0;
    node_type* pop_victim = &start;
    pop_victim = pop_victim->next[bottom_level].load();
    int x = pop_victim->node_value;
    while (pop_victim != &limit && !erase(x))
    {
      pop_victim = pop_victim->next[bottom_level].load();
      x = pop_victim->node_value;
    }
    if (pop_victim == &limit)
      return NULL;
    else
    {
      return pop_victim;
    }
  }
#endif


  template <class _Tp, class _Compare, class _Alloc, size_t _MAXLVL>
  bool skiplist<_Tp, _Compare, _Alloc, _MAXLVL> :: contains(const value_type& x)
  {
    // \pir \todo add memory management
    int        bottom_level = 0;
    bool       marked = false;
    node_type* pred = &start;
    node_type* curr = nullptr;
    node_type* succ = nullptr;

    for(int level = MAXLEVEL-1; level >= bottom_level; --level)
    {
      curr = pred->next[level].load();

      while (curr != &limit)
      {
        succ = curr->next[level].load();
        marked = curr->next[level].mark(std::memory_order_relaxed);

        while (marked)
        {
          curr = pred->next[level].load();

          if (curr == &limit)
            break;

          succ = curr->next[level].load();
          marked = curr->next[level].mark(std::memory_order_relaxed);
        }

        if (!(curr->node_value < x && curr != &limit)) break;

        pred = curr;
        curr = succ;
      }
    }

    return (curr->node_value == x);
  }

  //
  // quiescent methods

  template <class _Tp, class _Compare, class _Alloc, size_t _MAXLVL>
  typename skiplist<_Tp, _Compare, _Alloc, _MAXLVL>::iterator
  skiplist<_Tp, _Compare, _Alloc, _MAXLVL> :: qbegin()
  {
    // the first element is the one after start
    return iterator(start.next[0].load());
  }

  template <class _Tp, class _Compare, class _Alloc, size_t _MAXLVL>
  typename skiplist<_Tp, _Compare, _Alloc, _MAXLVL>::iterator
  skiplist<_Tp, _Compare, _Alloc, _MAXLVL> :: qend()
  {
    return iterator(&limit);
  }
} // namespace lockfree

#endif /* _SKIPLIST_HPP */
