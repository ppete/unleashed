/// \todo insert license

/// \file htm-skiplist.hpp
/// \author Created by Peter Pirkelbauer
///
/// \brief Loosely based on Herlihy and Shavit's fine grained locking skiplist implementation.
///        The Art of Multiprocessor Programming, 2012
///        The skiplist uses exclusively transactions to query and update the skiplist.
///        The implementation is lock-free; Any ongoing skiplist implementation
///        can be obstructed by a concurrently executing transaction.

#ifndef _HTMSKIPLIST_HPP

#define _HTMSKIPLIST_HPP 1

#include <memory>
#include <iterator>
#include <climits>
#include <vector>
#include <random>

#include "htm.hpp"
#include "htm-memory.hpp"
#include "bitutil.hpp"

namespace aux
{
  static const size_t MAGIC_INIT = 8171793;

  std::atomic<size_t> rndinit(7);

  thread_local std::ranlux48 rnd(rndinit.fetch_add(1) * MAGIC_INIT);

  static inline
  int bsr32(uint32_t s)
  {
    if (s == 0) return 0;

    static const int nobits = (sizeof(unsigned int) * CHAR_BIT);

    // note: use __builtin_clzl, and __builtin_clzll for long and long long
    return nobits - __builtin_clz(s);
  }

  constexpr
  size_t _next_power2(size_t N, size_t X)
  {
    return (X > 16) ? N : _next_power2(N | (N >> X), X << 1);
  }

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

  template <size_t M>
  static inline
  size_t log_rand()
  {
    const size_t X = next_power2(M);
    const size_t mask = (uint64_t(1) << (X)) - 1;
    const size_t r = (X - bsr32(rnd() & mask)) / 3;

    return (r < M) ? r : log_rand<M>();
  }
}

namespace htm
{
  template <class T>
  class SkipListNode
  {
    typedef SkipListNode<T>* NodePtr;

    public:
      typedef T value_type;

      NodePtr*              next;        ///< linkage levels; the number of
      const size_t          noLevels;    ///< number of linkage levels
                                         ///  computed randomly
      T                     nodeValue;   ///< user data

      // remove auto-generated ctors
      SkipListNode() = delete;
      SkipListNode(const SkipListNode&) = delete;
      SkipListNode& operator=(const SkipListNode&) = delete;

      /// Constructs a new skiplist node for an entry newentry; the node will
      /// participate on levelsLinked levels
      SkipListNode(const T& newentry, int levelsLinked)
      : next(new NodePtr[levelsLinked]), noLevels(levelsLinked), nodeValue(newentry)
      {
        assert(levelsLinked > 0);
      }

      NodePtr& operator[](size_t idx)
      {
        //~ htm_assert(next && idx < noLevels);

        if (!(next || idx < noLevels))
        {
          if (tx::transactional() == tx::active)
          {
            tx::end();
            std::cout << noLevels << " " << idx << " " << next << std::endl;
            assert(false);
          }
        }

        return next[idx];
      }

      ~SkipListNode()
      {
        delete[]next;
        //~ if (next)
        //~ {
          //~ NodePtr* old = next;
          //~ next = nullptr;
          //~ delete[]old;
        //~ }
      }

      /// number of linkage levels
      size_t levels()    const { return noLevels; }

      /// access to removed flag
      bool removed()     const { return next[0] == nullptr; }
      void mark_removed()      { next[0] = nullptr; }
  };

  template <class T>
  class SkipListIterator : std::iterator< std::forward_iterator_tag, T>
  {
    typedef std::iterator< std::forward_iterator_tag, T> base;

    public:
      typedef typename base::value_type        value_type;
      typedef typename base::reference         reference;
      typedef typename base::pointer           pointer;
      typedef typename base::difference_type   difference_type;
      typedef typename base::iterator_category iterator_category;

      explicit
      SkipListIterator(SkipListNode<T>* node)
      : skipNode(node)
      {
        assert(skipNode);
      }

      template <class U>
      friend
      bool operator==(const SkipListIterator<U>& lhs, const SkipListIterator<U>& rhs);

      template <class U>
      friend
      bool operator!=(const SkipListIterator<U>& lhs, const SkipListIterator<U>& rhs);

      reference operator*() const;
      SkipListIterator<T>& operator++();
      SkipListIterator<T>  operator++(int);

    private:
      SkipListNode<T>* skipNode;
  };

  template <class V>
  bool operator==(const SkipListIterator<V>& lhs, const SkipListIterator<V>& rhs)
  {
    return (lhs.skipNode == rhs.skipNode);
  }

  template <class V>
  bool operator!=(const SkipListIterator<V>& lhs, const SkipListIterator<V>& rhs)
  {
    return (lhs.skipNode != rhs.skipNode);
  }

  template <class V>
  typename SkipListIterator<V>::reference
  SkipListIterator<V>::operator*() const
  {
    return skipNode->nodeValue;
  }

  template <class T>
  SkipListIterator<T>& SkipListIterator<T> :: operator++()
  {
    skipNode = skipNode->next[0];

    assert(skipNode);
    return *this;
  }

  template <class T>
  SkipListIterator<T> SkipListIterator<T> :: operator++(int)
  {
    SkipListIterator<T> original(*this);
    skipNode = skipNode->next[0].load();
    assert(skipNode);
    return original;
  }

  template <int TXSIZE>
  static inline
  int _startofs()
  {
    return std::max(int(aux::rnd() % TXSIZE), 2);
  }

  // optimial startofs
  // std::atomic<int> ctr(0);

  template <class T>
  static inline
  int startofs(const counted<T>*)
  {
    // return 32 + (ctr.fetch_add(1, std::memory_order_relaxed) % 4) * 32;
    return tx::arch_tag::len;           /* fixed start offset */
    // return _startofs<tx::arch_tag::len>(); /* randomized start ofs */
  }

  static inline
  int startofs(const void*)
  {
    return tx::arch_tag::len;
  }

  template <class T>
  static inline
  size_t iterctr(const counted<T>*)
  {
    // return 32 + (ctr.fetch_add(1, std::memory_order_relaxed) % 4) * 32;
    // return tx::arch_tag::len;
    return 0;
  }

  static inline
  size_t iterctr(const void*)
  {
    return 1;
  }

  struct TXLength
  {
    static const size_t AVGSZ = 4;

    TXLength()
    : total(AVGSZ * tx::arch_tag::len), ctr(0)
    {
      recent[0] = recent[1] = recent[2] = recent[3] = tx::arch_tag::len;
    }

    size_t compute_next_len(size_t currlen, size_t iterctr)
    {
      size_t res = 0;

      if (iterctr > 0)
      {
        assert(currlen >= 2);

        total = total - recent[ctr] + currlen;
        recent[ctr] = currlen;
        ctr = (ctr + 1) % AVGSZ;
        res = total / AVGSZ;

        if (res == currlen)
        {
          res = std::min(2*tx::arch_tag::len, res + res/2);
        }
      }
      else
      {
        // ignore first iteration in ref_counter
        res = total / AVGSZ;
      }

      assert(res>=2);
      return res;
    }

    size_t  total;
    int64_t ctr;
    size_t  recent[4];
  };

  /// SkipList implementation
  /// \tparam _Tp element type
  /// \tparam _MAXLVL maximum number of levels in skiplist
  /// \tparam _Compare less-than operation
  /// \tparam _Alloc allocator; the default allocator uses a hazard pointer scheme
  template < class  _Tp,
             size_t _MAXLVL = 32,
             class  _Compare = std::less<_Tp>,
             class  _Alloc = pub_scan_manager<_Tp, std::allocator>
           >
  class SkipList
  {
    public:
      static const int    MAXLEVEL     = _MAXLVL;  // \todo set to reasonable amount
      static const int    NUMLOCALS    = 1;
      static const int    SZPINWALL    = MAXLEVEL + NUMLOCALS;

      static const size_t _MAXABORT    = 4;
      static const int    FIND_SUCCESS = -1;

      // \todo add boost concept checks
      typedef typename _Alloc::value_type                                   _Alloc_value_type;
      typedef std::allocator_traits<_Alloc>                                          _OrigAlloc_traits;
      typedef typename _OrigAlloc_traits::template rebind_alloc<_Tp>        _Tp_alloc_type;
      typedef typename _OrigAlloc_traits::template rebind_alloc<SkipListNode<_Tp> >  _Node_alloc_type;
      typedef std::allocator_traits<_Tp_alloc_type>                                  _Alloc_traits;

      typedef typename _Node_alloc_type::value_type                                  node_type;



      typedef typename _Node_alloc_type::pinguard                           PinGuard;

      struct find_data
      {
        node_type* pred;
        node_type* curr;
        int        level;
        int        levelFound;
        int        txlen;
#ifndef NDEBUG
        bool       faillate;
#endif /* NDEBUG */
      };

      typedef std::array<node_type*, MAXLEVEL>          level_array;

    public:
      typedef _Tp                                       value_type;
      typedef _Tp                                       key_type;

      typedef _Compare                                  value_compare;
      typedef _Compare                                  key_compare;

      typedef typename _Alloc_traits::pointer           pointer;
      typedef typename _Alloc_traits::const_pointer     const_pointer;
      typedef typename _Alloc_traits::value_type&       reference;
      typedef const typename _Alloc_traits::value_type& const_reference;

      typedef SkipListIterator<_Tp>                     iterator;
      typedef iterator                                  const_iterator; // \todo make real const iterator
      typedef typename iterator::difference_type        difference_type;

      typedef size_t                                    size_type;
      typedef _Alloc                                    allocator_type;

      explicit
      SkipList( const key_compare& comparator = key_compare(),
                const allocator_type& alloc = allocator_type()
              )
      : start(value_type(INT_MIN), _MAXLVL),
        limit(value_type(INT_MAX), _MAXLVL),
        lessThan(comparator),
        nodeAlloc(alloc),
        defaultLock(0)
      {
        static_assert(MAXLEVEL > 0, "MAXLEVEL > 0 required");

        for (int i = 0; i < MAXLEVEL; ++i)
        {
          start.next[i] = &limit;
          limit.next[i] = &limit; // link back to start?
        }
      }

      ~SkipList()
      {
        start.next = nullptr;
        limit.next = nullptr;
      }

      // \todo iterator find(const value_type& entry);

      /// \brief inserts a new entry into the list
      /// \return true, if the entry was inserted; false, otherwise.
      int insert(const value_type& newentry, size_t&);

      /// \brief removes entry from the list.
      /// \return true, if entry was in the list; false otherwise.
      int erase(const value_type& entry, size_t&);

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

      /// after an aborted transaction, we find a good starting point
      ///   for retrying.
      find_data
      find_restart_point(find_data data, level_array& preds);

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

    // private:
      /// Auxiliary Method to find the location of a given element
      inline
      find_data
      find(const value_type& x, level_array& preds, level_array& succs, find_data);

      /// Auxiliary Method to find the location of a given element
      find_data
      find_memsafe(const value_type& x, level_array& preds, level_array& succs, find_data);

      /// Auxiliary Method to find the location of a given element
      find_data find_init();

      node_type                start;
      node_type                limit;
      _Compare                 lessThan;
      _Node_alloc_type         nodeAlloc;
      std::atomic<size_t>      defaultLock;
      static thread_local      TXLength txestimater;
  };

  template <class _FindData>
  static inline
  void txlog(size_t abortcnt, _FindData data, std::string msg)
  {
    std::cerr << "abort: " << abortcnt
               << "   lv: " << std::setw(2) << data.level
               << "   cs: " << std::setw(2) << tx::abort_cause()
               << "   cd: " << std::setw(2) << tx::abort_code()
               << msg
               << std::endl;
  }

  template <class _FindData>
  static inline
  void txerr(size_t abortcnt, _FindData data, std::string msg, int exfailcode)
  {
    txlog(abortcnt, data, msg);

    throw exfailcode;
  }

  //
  // functions to model epoch_manager behavior that differs from other managers

  template <class _Tp, template <class> class _Alloc>
  static inline
  bool begin_tx_if_needed(epoch_manager<_Tp, _Alloc>* mgr)
  {
    // if the tx was canceled by another thread
    if (!mgr->active()) return false;

    return tx::begin();
  }

  template <class _Tp, template <class> class _Alloc>
  static inline
  void restart(epoch_manager<_Tp, _Alloc>* mgr)
  {
    // if the tx was canceled by another thread
    if (!mgr->active()) mgr->beginop(0);
  }

  template <class _Tp, template <class> class _Alloc>
  static inline
  void end_tx_if_needed(epoch_manager<_Tp, _Alloc>*)
  {
    htm_assert(tx::transactional() == tx::active);
    tx::end();
  }

  template <class _Tp, template <class> class _Alloc>
  void validate_transaction(epoch_manager<_Tp, _Alloc>* mgr)
  {
    mgr->validate_transaction();
  }

  static inline
  bool begin_tx_if_needed(void*)
  {
    return true;
  }

  static inline
  void end_tx_if_needed(void*) {}

  static inline
  void validate_transaction(void*) {}

  static inline
  void restart(void*) {}

#ifdef ACCELERATED_RESTART
static const bool accelerated_restart = true;
#else
static const bool accelerated_restart = false;
#endif

  template <class _Tp, size_t _MAXLVL, class _Compare, class _Alloc>
  typename SkipList<_Tp, _MAXLVL, _Compare, _Alloc> :: find_data
  SkipList<_Tp, _MAXLVL, _Compare, _Alloc> :: find_restart_point(find_data data, level_array& preds)
  {
    _Node_alloc_type alloc = getAllocator();

    if (accelerated_restart)
    {
      if (begin_tx_if_needed(&alloc))
      {
        while ((data.level < MAXLEVEL) && !preds[data.level]->removed()) ++data.level;
        end_tx_if_needed(&alloc);

        if (data.level < MAXLEVEL)
        {
          data.pred = data.curr = preds[data.level];
          return data;
        }
      }
    }

    restart(&alloc);
    return find_init();
  }

  template <class _Tp, size_t _MAXLVL, class _Compare, class _Alloc>
  int SkipList<_Tp, _MAXLVL, _Compare, _Alloc> :: insert(const value_type& x, size_t& abortctr)
  {
    const int maxLevel = aux::log_rand<MAXLEVEL>()+1;
    assert(maxLevel >= 1 && maxLevel <= MAXLEVEL);

    level_array      preds;
    level_array      succs;

    _Node_alloc_type alloc = getAllocator();
    size_t           myabort = 0;
    PinGuard         guard(alloc, SZPINWALL); // uses RAII to free alloc's pinwall
    find_data        data = find_init();

    for (;;)
    {
      int            txlen = startofs(&start);
      size_t         cntit = iterctr(&start);

      do
      {
        // store tx size
        data.txlen = txlen;

        if (tx::begin())
        {
          defaultLock.load(std::memory_order_relaxed);
          // invoke find transaction
          data = find_memsafe(x, preds, succs, data);
          tx::end();
          alloc.stats(txlen);
          txlen = txestimater.compute_next_len(txlen, cntit);
          ++cntit;
        }
        else
        {
          ++myabort;

          assert(data.pred && data.curr);
          assert(data.pred != &limit);

          // restart completely if remotely terminated
          if (!guard.active())
          {
            return insert(x, abortctr);
          }

          if (myabort >= _MAXABORT)
          {
            abortctr += myabort;
            myabort = 0;
            txlen = std::max(txlen / 2, 2);
          }
        }
      } while (data.level >= 0); // repeat while not yet found

      if (data.level == FIND_SUCCESS)
      {
        // return if key was found
        if (data.levelFound != -1) return -1;

        node_type* newNode = alloc.allocate(1);

        alloc.construct(newNode, x, maxLevel);

        // interlink new node
        for (int level = 0; level < maxLevel; ++level)
        {
          // newNode->next[level] = succs[level];
          (*newNode)[level] = succs.at(level);
        }

        htm_assert(tx::transactional() == tx::none);
        if (tx::begin())
        {
          defaultLock.load(std::memory_order_relaxed);
          for (int level = 0; level < maxLevel; ++level)
          {
            node_type* pred = preds.at(level);
            node_type* succ = succs.at(level);

            htm_assert(pred != &limit);

            if ( pred->removed() || (*pred)[level] != succ ) tx::abort<81>();
          }

          for (int level = 0; level < maxLevel; ++level)
          {
            node_type* pred = preds.at(level);

            pred->next[level] = newNode;
            // (*pred)[level] = newNode;
          }

          validate_transaction(&alloc);
          tx::end();
          assert(!data.faillate);
          abortctr += myabort;
          return maxLevel;
        }

        // restart completely if remotely terminated
        if (!guard.active())
        {
          return insert(x, abortctr);
        }

        alloc.deallocate(newNode, 1); // \todo consider introducing dealloc_nonshared
        abortctr += myabort + 1;
        myabort = 0;

        // if our anchor node is still part of the data structure
        //   we can retry from there.
        //   set level and pred accordingly
        data.level = maxLevel;
      } // end FIND_SUCCESS
      else
      {
        data.level += 2*MAXLEVEL;
      }

      data = find_restart_point(data, preds);
    }
  }

  template <class _Tp, size_t _MAXLVL, class _Compare, class _Alloc>
  typename SkipList<_Tp, _MAXLVL, _Compare, _Alloc>::find_data
  SkipList<_Tp, _MAXLVL, _Compare, _Alloc> :: find_init()
  {
    return find_data{ &start,                  /* pred */
                      &start,                  /* curr */
                      MAXLEVEL-1,              /* level */
                      -1,                      /* level_found */
                      tx::arch_tag::len        /* max tx size */
#ifndef NDEBUG
                      , false                  /* faillate */
#endif /* NDEBUG */
                    };
  }


  template <class _Tp, size_t _MAXLVL, class _Compare, class _Alloc>
  typename SkipList<_Tp, _MAXLVL, _Compare, _Alloc>::find_data
  SkipList<_Tp, _MAXLVL, _Compare, _Alloc> :: find(const value_type& x, level_array& preds, level_array& succs, find_data data)
  {
    htm_assert(data.curr);
    htm_assert(data.curr == data.pred);

    // make sure that the anchor point is still part of the data structure
    // -- if not communicate that through a negative value, we could use
    //    the level for better restart
    if (data.curr->removed()) { data.level -= 2*MAXLEVEL; return data; }

    // check consistency of the start node
    while (data.level >= 0)
    {
      while (data.curr != &limit && lessThan(data.curr->nodeValue, x))
      {
        data.pred = data.curr;

        // interrupt find when recommended transaction limit is reached.
        if (--data.txlen == 0)
        {
          preds[data.level] = data.pred;
          return data;
        }

        data.curr = static_cast<node_type*>(data.curr->next[data.level]);

        // inside a transaction, we should not see any removed nodes or nullptr
        htm_assert(data.curr && !data.curr->removed());
      }

      if (  data.levelFound == -1
         && data.curr != &limit
         && !lessThan(x, data.curr->nodeValue)
         )
      {
        data.levelFound = data.level;
      }

      preds.at(data.level) = data.pred;
      succs.at(data.level) = data.curr;
      --data.level;
      data.curr = data.pred;
    }

#if 0
    if        ((  data.pred == &start
              || data.pred == &limit
              || data.level <= 0
              || (*data.pred)[data.level] == nullptr
              || (*data.pred)[data.level-1] == nullptr
              || (*data.pred)[data.level]->nodeValue >= data.pred->next[data.level-1]->nodeValue
              ) == false)
    {
      data.faillate = true;
    }

    if        ((  data.curr == &start
              || data.curr == &limit
              || data.level <= 0
              || (*data.curr)[data.level] == nullptr
              || (*data.curr)[data.level-1] == nullptr
              || (*data.curr)[data.level]->nodeValue >= data.curr->next[data.level-1]->nodeValue
              ) == false)
    {
      data.faillate = true;
    }
#endif

    return data;
  }

  template <class _Tp, size_t _MAXLVL, class _Compare, class _Alloc>
  typename SkipList<_Tp, _MAXLVL, _Compare, _Alloc>::find_data
  SkipList<_Tp, _MAXLVL, _Compare, _Alloc> :: find_memsafe(const value_type& x, level_array& preds, level_array& succs, find_data data)
  {
    data = find(x, preds, succs, data);

    if (data.level >= FIND_SUCCESS)
    {
      // pin all predecessors that have been identified
      _Node_alloc_type alloc = getAllocator();
      int              lower_bound = data.level;

      if (lower_bound < 0) lower_bound = 0;
      alloc.txpin(preds.begin() + lower_bound, preds.end(), &start);
    }

    return data;
  }

  template <class _Tp, size_t _MAXLVL, class _Compare, class _Alloc>
  int SkipList<_Tp, _MAXLVL, _Compare, _Alloc> :: erase(const value_type& x, size_t& abortctr)
  {
    _Node_alloc_type alloc   = getAllocator();
    size_t           myabort = 0;
    PinGuard         guard(nodeAlloc, SZPINWALL);
    level_array      preds;
    level_array      succs;
    find_data        data    = find_init();

    for (;;)
    {
      int              txlen  = startofs(&start);
      int              cntit  = iterctr(&start);

      do
      {
        // store tx size
        data.txlen = txlen;

        if (tx::begin())
        {
          defaultLock.load(std::memory_order_relaxed);
          // invoke find transaction
          data = find_memsafe(x, preds, succs, data);
          tx::end();
          alloc.stats(txlen);
          txlen = txestimater.compute_next_len(txlen, cntit);
          ++cntit;
        }
        else
        {
          ++myabort;

          if (!guard.active())
          {
            return erase(x, abortctr);
          }

          if (myabort >= _MAXABORT)
          {
            abortctr += myabort;
            myabort = 0;
            txlen = std::max(txlen / 2, 2);
          }
        }
      } while (data.level >= 0); // repeat while not yet found

      if (data.level == FIND_SUCCESS)
      {
        abortctr += myabort;
        myabort = 0;
        // if node is not found in the list
        if (data.levelFound == -1) return -1;

        node_type* victim          = succs.at(data.levelFound);
        assert(victim);

        assert(tx::transactional() == tx::none);
        if (tx::begin())
        {
          // in case somebody else is removing this node
          if (victim->removed())
          {
            tx::end();
            return -1;
          }

          defaultLock.load(std::memory_order_relaxed);
          int maxLevel        = victim->levels()-1;

          for (int level = 0; level <= maxLevel; ++level)
          {
            node_type* pred = preds.at(level);
            htm_assert(pred);

            // validate pred
            // if (pred->next[level] != victim) tx::abort<79>();
            if ((*pred)[level] != victim) tx::abort<79>();
            if (pred->removed()) tx::abort<78>();

            node_type* succ = static_cast<node_type*>((*victim)[level]);

            htm_assert(succ);
            if (succ->removed()) tx::abort<80>();
          }

          for (int level = maxLevel; level >= 0; --level)
          {
            node_type* pred = preds.at(level);
            node_type* succ = static_cast<node_type*>((*victim)[level]);

            (*pred)[level] = succ;
          }

          // kill transactions that rely on victim being current in the data-structure
          victim->mark_removed();
          validate_transaction(&alloc);
          tx::end();

          // validate that we have seen a properly interlinked victim
          assert(data.levelFound == maxLevel);
          assert(!data.faillate);

          // deallocate node
          alloc.deallocate(victim, 1);
          abortctr += myabort;
          return maxLevel+1;
        }

        ++myabort;
        data.level = data.levelFound;

        if (!guard.active())
        {
          return erase(x, abortctr);
        }
      }
      else
      {
        ++myabort;
        data.level += 2*MAXLEVEL;
      }

      data = find_restart_point(data, preds);
    }

    assert(false);
    return -1;
  }




  template <class _Tp, size_t _MAXLVL, class _Compare, class _Alloc>
  bool SkipList<_Tp, _MAXLVL, _Compare, _Alloc> :: contains(const value_type& x)
  {
    std::array<node_type*, MAXLEVEL> preds;
    std::array<node_type*, MAXLEVEL> succs;
    PinGuard                         guard(nodeAlloc, SZPINWALL); // uses RAII to free alloc's pinwall
    const int                        lFound = find(x, preds, succs);

    // \mo
    // sequentially consistent, since we do not acquire any locks
    return (  lFound != -1
           && succs[lFound]->fullyLinked.load()
           && !succs[lFound]->marked.load()
           );
  }


  // \pp \todo not sure if needed; unsafe for concurrent environment
  template <class _Tp, size_t _MAXLVL, class _Compare, class _Alloc>
  size_t SkipList<_Tp, _MAXLVL, _Compare, _Alloc> :: qsize() // \todo const
  {
    iterator aa = qbegin();
    iterator zz = qend();
    size_t   cnt = 0;

    while (aa != zz) { ++cnt; ++aa; }

    return cnt;
  }

  template <class _Tp, size_t _MAXLVL, class _Compare, class _Alloc>
  SkipListIterator<_Tp> SkipList<_Tp, _MAXLVL, _Compare, _Alloc> :: qbegin()
  {
    // \pp needed?
    // wait for it to not be marked
    return iterator(start.next[0]);
  }

  template <class _Tp, size_t _MAXLVL, class _Compare, class _Alloc>
  SkipListIterator<_Tp> SkipList<_Tp, _MAXLVL, _Compare, _Alloc> :: qend()
  {
    return iterator(&limit);
  }


  template <class _SkipNode>
  std::vector<_SkipNode*> create_links(_SkipNode& start)
  {
    typedef std::vector<_SkipNode*> node_vector;
    // return node_vector( &start.next[0], &start.next[start.levels()] );

    node_vector res;

    for (size_t i = 0; i < start.levels(); ++i)
      res.push_back(static_cast<_SkipNode*>(start.next[i]));

    return res;
  }


  template <class _SkipNode, class LessThan>
  bool check_links(std::vector<_SkipNode*> links, _SkipNode& lhs, _SkipNode& rhs, LessThan lessThan, size_t hops)
  {
    bool valid = lessThan(lhs.nodeValue, rhs.nodeValue);

    if (lhs.removed())
    {
      std::cout << "removed " << lhs.nodeValue << std::endl;
      valid = false;
    }

    for (size_t i = 0; i < rhs.levels(); ++i)
    {
      if (  (links[i] != &rhs)
         || (links[i]->levels() <= i)
         )
      {
        std::cout << "i : " << i
                  << "   |l.lv| = " << lhs.levels() << (lhs.levels() <= i ? "!" : "" ) << " v = " << lhs.nodeValue
                  << "   |r.lv| = " << rhs.levels() << (rhs.levels() <= i ? "!" : "" ) << " v = " << rhs.nodeValue
                  << "   |off!| = " << links[i]->levels() << (links[i]->levels() <= i ? "!" : "" ) << " v = " << links[i]->nodeValue
                  << "   hops = " << hops
                  << std::endl;

        valid = false;
      }
    }

    return valid;
  }

  template <class _SkipNode>
  std::vector<_SkipNode*>
  update_links(std::vector<_SkipNode*> links, _SkipNode& n)
  {
    for (size_t i = 0; i < n.levels(); ++i)
    {
      links[i] = static_cast<_SkipNode*>(n.next[i]);
    }

    return links;
  }

  template <class _SkipNode, class LessThan>
  bool skiplist_check(_SkipNode* aa, _SkipNode* zz, LessThan lessThan)
  {
    std::vector<_SkipNode*> links = create_links(*aa);

    _SkipNode* prev = aa;
    _SkipNode* curr = static_cast<_SkipNode*>(aa->next[0]);
    size_t     hops = 0;
    bool       vald = true;

    while (prev != zz)
    {
      vald = check_links(links, *prev, *curr, lessThan, ++hops) && vald;
      links = update_links(links, *curr);

      prev = curr;
      curr = static_cast<_SkipNode*>(curr->next[0]);
    }

    return vald;
  }

  template <class _SkipNode, class LessThan>
  bool skiplist_check(_SkipNode& aa, _SkipNode& zz, LessThan lessThan)
  {
    return skiplist_check(&aa, &zz, lessThan);
  }

  template <class _Tp, size_t _MAXLVL, class _Compare, class _Alloc>
  thread_local
  TXLength SkipList<_Tp, _MAXLVL, _Compare, _Alloc> :: txestimater;
}

#endif /* _HTMSKIPLIST_HPP */




#if 0
      if (myabort >= _MAXABORT - 5)
      {
        std::cerr << "----- ----- -----"
                  << "abort: " << myabort
                  << "   lv: " << std::setw(2) << data.level
                  << "   cs: " << std::setw(2) << tx::abort_cause()
                  << "   cd: " << std::setw(2) << tx::abort_code()
                  << "   lf: " << (data.levelFound)
                  << "   |v|: " << (succs[data.levelFound]->levels() - 1)
                  << "   e"
                  << std::endl;

        std::cerr << "v = " << succs[data.levelFound]
                  << "   (" << succs[data.levelFound]->nodeValue << ")"
                  << ( succs[data.levelFound]->removed ? "DEL" : "" )
                  << std::endl ;

        for (int level = (succs[data.levelFound]->levels() - 1); level >= 0; --level)
        {
          node_type* pred = preds[level];

          // validate pred
          std::cerr << "p = " << pred
                    << "   (" << pred->nodeValue << ")"
                    << ( pred->removed ? "DEL" : "" )
                    << std::endl;

          std::cerr << "p-> " << pred->next[level]
                    << "   (" << pred->next[level]->nodeValue << ")"
                    << ( pred->next[level]->removed ? "DEL" : "" )
                    << "   @" << level
                    << std::endl;

          std::cerr << "v-> " << victim->next[level]
                    << "   (" << victim->next[level]->nodeValue << ")"
                    << ( victim->next[level]->removed ? "DEL" : "" )
                    << "   @" << level
                    << std::endl;
        }
      }
#endif
