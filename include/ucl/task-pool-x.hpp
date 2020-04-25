#ifndef _TASK_POOL_X_HPP
#define _TASK_POOL_X_HPP 1

#include "ucl/pmemory.hpp"
#include "ucl/task.hpp"
#include "ucl/task-simple-queue.hpp"

namespace ucl
{
  namespace
  {
    inline
    bool is_valid_pointer(const void* ptr)
    {
      return (  (ptr != nullptr)
             && ( (reinterpret_cast<intptr_t>(ptr) & 1) == 0 )
             );
    }
  }
  
  /// \private
  /// \brief Data block in a doubly linked queue of FIFO queues.
  ///        doubleQ supports a single producer, multiple consumers
  ///        tail->block<->block<->block<->head
  ///        with the producer accessing the tail and thieves accessing
  ///        the head.
  ///
  /// \details A data block holds BLKSZ elements at a time.
  ///          The data block is NOT cyclic, hence we get away with
  ///          relaxed memory operations. It is the task of the memory
  ///          management to guarantee that the deq operation completes
  ///          before the block can be recycled.
  /// \tparam T the task being stored
  template <class T>
  struct doubleQ
  {
      /// type for head and tail position within block
      typedef ucl::aligned_atomic_type<nat_t, CACHELINESZ>      counter_t;
  
      /// pointer to next data block
      typedef ucl::aligned_atomic_type<doubleQ<T>*, CACHELINESZ> next_t;
  
      /// task type
      typedef T                                                  value_t;
  
      static const nat_t BLK = ucl::DEFAULT_BLOCK_SIZE; ///< number of elements in a block
  
      counter_t   hd;         ///< block internal head
      counter_t   tl;         ///< block internal tail
      next_t      next;       ///< next block
      doubleQ<T>* prev;       ///< only accessed by owner
  
      /// elements in the block
      /// \note in order to avoid interference from ctor and dtor calls
      ///       by classes with non-trivial constructor and deconstructor, the
      ///       block manages construction and deconstruction of its task objects.
      ///       Thus, it just reserves storage for BLK elements, and calls ctor
      ///       and dtor upon enqueue and dequeue (move based when possible).
      char       data[BLK*sizeof(value_t)];
  
      /// constructs data block with 0 elements
      doubleQ()
      : hd(0), tl(0), next(nullptr), prev(nullptr)
      {}
  
      /// constructs a data block with 1 element
      /// \param el the first task
      explicit
      doubleQ(value_t&& el)
      : hd(0), tl(1), next(nullptr), prev(nullptr)
      {
        new (at(0)) value_t (std::move(el));
      }
  
      /// constructs a data block with 1 element
      /// \param el the first task
      explicit
      doubleQ(const value_t& el)
      : hd(0), tl(1), next(nullptr), prev(nullptr)
      {
        new (at(0)) value_t (el);
      }
      
      /// tests if this block is empty
      bool is_empty() const
      {
        return hd.val.load(std::memory_order_relaxed) >= BLK;
      }
  
      /// computes the address of element idx
      value_t* at(nat_t idx)
      {
        char* res = data + idx*sizeof(T);
  
        assert(res >= data && res <= data + (BLK-1)*sizeof(T));
        return reinterpret_cast<value_t*>(res);
      }
  
      /// tests if the data block can store more tasks
      bool has_space()
      {
        // \mo relaxed, b/c tl is only modified by this thread
        return tl.val.load(std::memory_order_relaxed) < BLK;
      }
  
      /// copy-enqueues a new element in this block
      /// \pre hasSpace() == true
      /// \pre only executed by the queue owner
      void enq(const value_t& el)
      {
        // \mo relaxed, b/c tl is only modified by this thread
        const nat_t tail = tl.val.load(std::memory_order_relaxed);
  
        assert(tail < BLK);
        new (at(tail)) value_t (el); // in-place constructions
  
        // \mo release b/c data[tail] is published
        tl.val.store(tail+1, std::memory_order_release);
      }
  
      /// move-enqueues a new element in this block
      /// \pre hasSpace() == true
      /// \pre only executed by the queue owner
      void enq(value_t&& el)
      {
        // \mo relaxed, b/c tl is only modified by this thread
        const nat_t tail = tl.val.load(std::memory_order_relaxed);
  
        assert(tail < BLK);
        new (at(tail)) value_t (std::move(el)); // in-place constructions
  
        // \mo release b/c data[tail] is published
        tl.val.store(tail+1, std::memory_order_release);
      }
  
      static
      doubleQ<T>* deq(doubleQ<T>* self, void* rawmem);
  
      /// dequeues a task and moves it into res
      /// \result the block from where the element was dequeued;
      ///         nullptr if unsuccessful.
      doubleQ<T>* deq(void* rawmem)
      {
        return deq(this, rawmem);
      }
  
      template <class Alloc>
      static
      T* deq_steal(doubleQ<T>* self, void* rawmem, const int_fast8_t tries, Alloc alloc);
  
      /// steals  a task and moves it into res
      /// \param  res the place where the task will be stored
      /// \param  tries number of attempts. Used to manage lock-free contention
      ///         by giving up after tries unsuccessful attempts to steal a task.
      /// \param  alloc the memory manager in use; needed to secure the next block
      ///         if traversal is needed.
      /// \result the block from where the element was dequeued;
      ///         nullptr if unsuccessful.
      /// \pre    this is pinned, if a finegrain (i.e., publishAndScan) is used
      /// \post   this is unpinned
      /// \note   publishAndScan is supported, though finegrain memory management
      ///         does not bring anything to the table for this task implementation.
      template <class Alloc>
      T* deq_steal(void* rawmem, const int_fast8_t tries, Alloc alloc)
      {
        return deq_steal(this, rawmem, tries, alloc);
      }
      
    private:
      doubleQ(const doubleQ&)            = delete;
      doubleQ(doubleQ&&)                 = delete;
      doubleQ& operator=(const doubleQ&) = delete;
      doubleQ& operator=(doubleQ&&)      = delete;
  };

  template <class T>
  doubleQ<T>*
  doubleQ<T>::deq(doubleQ<T>* self, void* rawmem)
  {
    do
    {
      // \mo relaxed, b/c hd does not publish anything
      nat_t head = self->hd.val.load(std::memory_order_relaxed);
      // \mo relaxed, b/c we are on the same thread
      nat_t tail = self->tl.val.load(std::memory_order_relaxed);

      assert(tail <= BLK);
      while (head < tail)
      {
        assert(head < BLK);

        // \mo relaxed (succ) since the data is never rewritten
        // \mo relaxed (fail) b/c nothing happened
        if (self->hd.val.compare_exchange_strong(head, head+1, std::memory_order_relaxed, std::memory_order_relaxed))
        {
          assert(head+1 <= tail);
          value_t* obj = self->at(head);

          // mv and deconstruct obj
          new (rawmem) value_t(std::move(*obj));
          obj->~value_t();
          return self;
        }
      }

      self = self->prev;
    } while (self);

    return nullptr;
  }

  template <class T>
  template <class Alloc>
  T*
  doubleQ<T>::deq_steal(doubleQ<T>* self, void* rawmem, const int_fast8_t tries, Alloc alloc)
  {
    assert(tries > 0);
    assert(is_valid_pointer(self));

    do
    {
      assert(is_valid_pointer(self));
      
      int_fast8_t attempts = tries;
      // \mo relaxed, b/c hd does not publish anything
      const nat_t hd2  = self->hd.val.load(std::memory_order_relaxed);
      nat_t       head = hd2;
      // \mo acquire, b/c we read from data[i], where i < tl
      nat_t       tail = self->tl.val.load(std::memory_order_acquire);

      assert(head <= BLK);
      
      // give up under high contention
      while (head != BLK)
      {
        assert (head < BLK);
        if (head >= tail) return nullptr;

        // \mo relaxed (succ) since the data is never rewritten
        // \mo relaxed (fail) b/c nothing happened
        if (self->hd.val.compare_exchange_strong(head, head+1, std::memory_order_relaxed, std::memory_order_relaxed))
        {
          assert(head < tail);
          value_t* obj = self->at(head);

          new (rawmem) value_t (std::move(*obj));
          obj->~value_t();
          return static_cast<T*>(rawmem);
        }

        if (--attempts == 0) return nullptr;
      }

      doubleQ<T>* tmp = self;

      // \mo consume, b/c we want to see initialized memory
      self = alloc.template pin<std::memory_order_consume>(self->next.val);

      alloc.unpin(tmp, -1);
    } while (self);

    return nullptr;
  }


  /// \private
  /// Main FIFO queue, built upon doubleQ blocks + status data for threads
  /// \tparam T task type
  /// \tparam Alloc the memory manager as defined in pmemory.hpp
  /// \details
  ///    Each thread has its own dequeQ object. The dequeQ object also holds
  ///    some status information for its owning thread.
  template <class T, class UclAlloc>
  struct dequeQ
  {
      /// Memory manager type
      typedef UclAlloc                                                       orig_alloc_type;
      typedef std::allocator_traits<UclAlloc>                                orig_alloc_traits;
      typedef typename orig_alloc_traits::template rebind_alloc<doubleQ<T> > allocator_type;
      typedef T                                                              value_type;

      /// Number of concurrently pinned blocks (this and next)
      static const nat_t                                 SZPINWALL = 2;

    //~ private:
      /// The allocator/memory manager
      allocator_type                                     nodealloc;

      /// tail (enqueue-side) of the queue; only accessed by owner
      doubleQ<T>*                                        tail; // only accessed by owner

      /// head (dequeue-side) of the queue; concurrently accessed by work stealers.
      ucl::aligned_atomic_type<doubleQ<T>*, CACHELINESZ> head;

      /// number of tasks stolen from this thread
      nat_t                                              stolen;

      /// number of tasks enqueued by this thread
      nat_t                                              tasks_made;

      /// number of completed tasks (incl. tasks solved by other threads, but
      ///   not tasks that this thread stole).
      nat_t                                              tasks_done;

      /// number of stolen and completed tasks; updated by the work stealers.
      ucl::aligned_atomic_type<nat_t, CACHELINESZ>       steals;

      /// estimate of available tasks
      ucl::aligned_atomic_type<nat_t, CACHELINESZ>       numtasks;

    public:
      /// default ctor initializes an empty Q
      explicit
      dequeQ(orig_alloc_type alloc = orig_alloc_type())
      : nodealloc(alloc), tail(nullptr), head(nullptr), stolen(0),
        tasks_made(0), tasks_done(0), steals(0), numtasks(0)
      {
        nodealloc.initialize_if_needed();

        tail = nodealloc.allocate(1);
        new (tail)doubleQ<T>();
        head.val.store(tail);
      }

      /// retrieve allocator
      allocator_type get_allocator()
      {
        return nodealloc;
      }

      /// returns an estimated number of available tasks
      nat_t estimate_size() const
      {
        return numtasks.val.load(std::memory_order_relaxed);
      }

      /// allocates a new doubleQ block
      doubleQ<T>* new_block(T&& el)
      {
        allocator_type alloc = get_allocator();
        doubleQ<T>*    tmp   = alloc.allocate(1);

        new (tmp)doubleQ<T>(std::move(el));
        return tmp;
      }

      /// copy-enqueues a new element
      void enq(const T& el)
      {
        ++tasks_made;

        if (tail->has_space())
        {
          tail->enq(el);
          return;
        }

        doubleQ<T>* tmp = new_block(std::move(T(el)));

        tmp->prev = tail;
        numtasks.val.store(tasks_made - tasks_done, std::memory_order_relaxed);

        // \mo release, b/c new_block was initialized
        tail->next.val.store(tmp, std::memory_order_release);
        tail = tmp;
      }

      /// move-enqueues a new element
      void enq(T&& el)
      {
        ++tasks_made;

        if (tail->has_space())
        {
          tail->enq(std::move(el));
          return;
        }

        doubleQ<T>* tmp = new_block(std::move(el));

        numtasks.val.store(tasks_made - tasks_done, std::memory_order_relaxed);

        // \mo release, b/c new_block was initialized
        tail->next.val.store(tmp, std::memory_order_release);
        tmp->prev = tail;
        tail = tmp;
      }

      /// releases full tiles at head and tail
      void release_mem(doubleQ<T>* limit);
      
      /// releases memory to memory manager w/o actually freeing it 
      void release_mem_core(doubleQ<T>* limit);

      /// dequeues an element and stores result into res
      /// \param  res  where the result will be stored
      /// \result true if successful, false otherwise
      T* deq(void* rawmem)
      {
        assert(tail);

        // attempt to dequeue from own queue
        doubleQ<T>* actl = tail->deq(rawmem);

        // unsuccessful?
        if (actl == nullptr) return nullptr;

        // all blocks that have become empty can be freed
        if (tail != actl)
        {
          release_mem(actl);

          // update available task counter
          numtasks.val.store(tasks_made - tasks_done, std::memory_order_relaxed);
        }

        return static_cast<T*>(rawmem);
      }

      /// attempts to steal a task from some other thread's queue
      /// \param  res   where the result will be stored
      /// \param  tries number of attempts before yielding to contention
      /// \result true if successful, false otherwise
      T* deq_steal(void* rawmem, int_fast8_t tries)
      {
        assert(tries);

        allocator_type alloc = get_allocator();

        // \mo consume, b/c we need to see memory initialized
        doubleQ<T>*    curr = alloc.template pin<std::memory_order_consume>(head.val);
        assert(is_valid_pointer(curr));

        return curr->deq_steal(rawmem, tries, alloc);
      }

      /// records finished task
      void task_completed() { ++tasks_done; }

      /// records finished stolen task
      void stolen_task_completed()
      {
        steals.val.fetch_add(1, std::memory_order_relaxed);
      }

      /// does any created task (stolen or not) remain unfinished?
      bool has_unfinished_tasks()
      {
        nat_t prev_stolen = stolen;

        stolen = steals.val.load(std::memory_order_relaxed);

        assert(tasks_made - tasks_done >= stolen - prev_stolen);
        tasks_done += (stolen - prev_stolen);

        return tasks_made > tasks_done;
      }

      /// quiescent release of memory
      /// \details
      ///    DO NOT CALL when concurrent operations are ongoing!
      void qshutdown(bool isMainThread)
      {
        // release head will not work, since tail might not be empty yet
        release_mem_core(tail);
        
        get_allocator().deallocate(tail, 1);
        get_allocator().qrelease_memory();
        get_allocator().qdestruct(isMainThread);
      }

    private:
      dequeQ(const dequeQ&)            = delete;
      dequeQ(dequeQ&&)                 = delete;
      dequeQ& operator=(const dequeQ&) = delete;
      dequeQ& operator=(dequeQ&&)      = delete;
  };
  
  template <class T, class UclAlloc>
  void dequeQ<T, UclAlloc>::release_mem_core(doubleQ<T>* limit)
  { 
    assert(limit != nullptr);    
    allocator_type alloc = get_allocator();

    while (tail != limit && tail->is_empty())
    {
      doubleQ<T>* temp = tail;

      tail = tail->prev;
      tail->next.val.store(nullptr, std::memory_order_relaxed);
      
      // do tmp's links have to be reset?
      temp->prev = nullptr;
      temp->next.val.store(nullptr, std::memory_order_relaxed);

      alloc.deallocate(temp, 1);
    }

    // \mo intra thread
    doubleQ<T>* curr      = head.val.load(std::memory_order_relaxed);
    std::vector<doubleQ<T>*> releaseable;
    
    while (curr != limit && curr->is_empty())
    {
      releaseable.push_back(curr);
      curr = curr->next.val.load(std::memory_order_relaxed);
    }
    
    if (releaseable.size())
    {    
      // \mo store head and then prevent reordering with following deallocate calls
      head.val.store(curr, std::memory_order_relaxed);
      std::atomic_thread_fence(std::memory_order_acq_rel);
      
      curr->prev = nullptr;
      for (doubleQ<T>* elem : releaseable)
      {
        alloc.deallocate(elem, 1);
      }
    }
  }
  
  template <class T, class UclAlloc>
  void dequeQ<T, UclAlloc>::release_mem(doubleQ<T>* limit)
  {
    release_mem_core(limit);
    
    get_allocator().release_memory_if_needed();

#if INFORMED_STEALING
    // update completed stolen tasks
    has_unfinished_tasks();
#endif /* INFORMED_STEALING */
  }

  template <class T>
  using dbg_task_allocator = lockfree::epoch_manager<T, std::allocator>;
  //~ using dbg_task_allocator = lockfree::just_alloc<T, std::allocator>;

  template <class T>
  using deque_like_pool = ucl::perthread_pool<ucl::dequeQ<T, dbg_task_allocator<T> > >;
  //~ using deque_like_pool = ucl::perthread_pool<ucl::dequeQ<T, default_task_allocator<T> > >;


  template <class T>
  using default_pool = deque_like_pool<T>;

  //~ template <class T>
  //~ using default_pool = fifo_queue_pool<T>;
  //template <class T>
  //using default_pool = single_q_pool<locking::queue<T, ucl::ttas_lock_backoff_default> >;

  //~ template <class T>
  //~ using default_pool = single_q_pool<locking::queue<T> >;

  //template <class T>
  //using default_pool = single_q_pool<locking::queue<T, ucl::ttas_lock, ucl::lock_elision_guard> >;

  //~ template <class T>
  //~ using default_pool
    //~ = single_q_pool<lockfree::queue<T, lockfree::epoch_manager<T, std::allocator> > >;

  //~ template <class T>
  //~ using default_pool
    //~ = single_q_pool<lockfree::queue<T, lockfree::pub_scan_manager<T, std::allocator> > >;


  template <class F, class T>
  static inline
  auto execute_tasks_x(nat_t numthreads, F fun, T&& task)
       -> decltype(execute_pool(fun, *new default_pool<T>(numthreads, std::move(task))))
  {
    default_pool<T> pool(numthreads, std::move(task));

    return execute_pool(fun, pool);
  }
}

#endif /* _TASK_POOL_X_HPP */
