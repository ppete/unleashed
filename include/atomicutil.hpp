#ifndef _ATOMICUTIL_HPP
#define _ATOMICUTIL_HPP

/// File providing various atomic utility functions and classes
/// \author Peter Pirkelbauer
/// \email  pirkelbauer@uab.edu

#include <atomic>
#include <cassert>

#ifndef CACHELINESZ
  #if defined(_ARCH_PPC64)
    #define CACHELINESZ 128
  #else
    #define CACHELINESZ 64
  #endif
#endif

namespace uab
{
  template <class T, size_t ALIGNSZ>
  struct aligned_type
  {
    static_assert(ALIGNSZ >= sizeof(T), "alignment/size mismatch");
    
    T              val;
    char           x[ALIGNSZ-sizeof(T)];

    aligned_type(T v)
    : val(v)
    {}

    aligned_type()
    : val(T())
    {}
  };

  template <class T, size_t ALIGNSZ>
  struct aligned_atomic_type
  {
    static_assert(ALIGNSZ >= sizeof(std::atomic<T>), "alignment/size mismatch");

    std::atomic<T> val;
    char           x[ALIGNSZ-sizeof(std::atomic<T>)];

    aligned_atomic_type(T v)
    : val(v)
    {}

    aligned_atomic_type()
    : val(T())
    {}
  };


  /// Class providing similar functionality to Java's AtomicMarkableReference
  /// \tparam T underlying type; class represents markable T*
  /// \tparam MARKABLE_BITS number of expected markable bits
  template <class T, size_t MARKABLE_BITS = 1>
  struct MarkablePointer
  {
      // \todo replace uintptr_t w/ type mapped from T* to make sure that
      //       the sizes are compatible
      typedef uintptr_t                numptr_type; ///< unsigned int type that can hold T*
      typedef uintptr_t                mark_type;   ///< unsigned int type that can hold the marked bits
      typedef std::pair<T*, mark_type> state_type;  ///< pair of pointer-value + mark bits

    private:
      static const size_t MAX_MARKED = 1 << (MARKABLE_BITS+1);
      static const size_t MARKBITS   = MAX_MARKED-1;

      /// internal function to convert from T* to the internal
      ///   representation numptr_type
      static
      numptr_type toNum(T* ptr)
      {
        return reinterpret_cast<numptr_type>(ptr);
      }

      /// internal function to convert from the internal representation numptr_type
      ///   to T*
      static
      T* toPtr(numptr_type val)
      {
        return reinterpret_cast<T*>(val);
      }

    public:
      /// Constructs a new MarkablePointer with provided values.
      /// \details Constructs an unmarked nullptr by default.
      explicit
      MarkablePointer(T* ptr = nullptr, mark_type markflags = 0)
      : val(0)
      {
        assert((markflags < MAX_MARKED) && ((toNum(ptr) & MARKBITS) == 0));

        val.store(toNum(ptr) | markflags); // \mo
      }

      /// Constructs a new MarkablePointer from a markable pointer state.
      explicit
      MarkablePointer(state_type s)
      : MarkablePointer(s.first, s.second)
      {}

      /// Reads the current pointer.
      /// \note comparable to AtomicMarkableReference.getReference
      T* load(std::memory_order order = std::memory_order_seq_cst) const
      {
        return toPtr(val.load(order) & ~MARKBITS);
      }

      /// Sets the current value and markbits.
      /// \note comparable to AtomicMarkableReference.set
      void store(T* ptr, mark_type markflags = 0, std::memory_order order = std::memory_order_seq_cst)
      {
        assert((markflags < MAX_MARKED) && ((toNum(ptr) & MARKBITS) == 0));

        return val.store(toNum(ptr) | markflags, order);
      }

      /// Attempts to mark the reference.
      /// \details Replaces the current mark of this reference with the new
      ///          mark that is provided using compare_exchange_weak.
      /// \return  true, iff the reference's mark was successfully updated.
      /// \pre     expectedUnmarked is an unmarked pointer
      /// \note    comparable to AtomicMarkableReference.attemptMark
      bool mark_weak( T* expectedUnmarked,
                      mark_type markflags = 1,
                      std::memory_order succ = std::memory_order_seq_cst,
                      std::memory_order fail = std::memory_order_seq_cst
                    )
      {
        assert(  (markflags < MAX_MARKED)
              && ((toNum(expectedUnmarked) & MARKBITS) == 0)
              );

        numptr_type ov = val.load(std::memory_order_relaxed);
        numptr_type nv = toNum(expectedUnmarked) | markflags;

        return val.compare_exchange_weak(ov, nv, succ, fail);
      }

      /// Attempts to mark the reference.
      /// \details Replaces the current mark of this reference with the new
      ///          mark that is provided using compare_exchange_strong.
      /// \return  true, iff the reference's mark was successfully updated.
      /// \pre     expectedUnmarked is an unmarked pointer
      bool mark_strong( T* expectedUnmarked,
                        mark_type markflags = 1,
                        std::memory_order succ = std::memory_order_seq_cst,
                        std::memory_order fail = std::memory_order_seq_cst
                      )
      {
        assert(  (markflags < MAX_MARKED)
              && ((toNum(expectedUnmarked) & MARKBITS) == 0)
              );

        numptr_type ov = val.load(std::memory_order_relaxed);
        numptr_type nv = toNum(expectedUnmarked) | markflags;

        return val.compare_exchange_strong(ov, nv, succ, fail);
      }

      /// Weak compare_and_set based on old and new pairs of pointer and mark values.
      /// \note Comparable to AtomicMarkableReference.WeakCompareAndSet
      bool compare_exchange_weak( T*&               oldptr,
                                  T*                newptr,
                                  mark_type&        oldmark,
                                  mark_type         newmark,
                                  std::memory_order succ = std::memory_order_seq_cst,
                                  std::memory_order fail = std::memory_order_seq_cst
                                )
      {
        assert(  (oldmark < MAX_MARKED)
              && ((toNum(oldptr) & MARKBITS) == 0)
              && (newmark < MAX_MARKED)
              && ((toNum(newptr) & MARKBITS) == 0)
              );

        numptr_type       oldnum = reinterpret_cast<numptr_type>(oldptr) | oldmark;
        const numptr_type newnum = reinterpret_cast<numptr_type>(newptr) | newmark;
        const bool     res    =  val.compare_exchange_weak(oldnum, newnum, succ, fail);

        oldptr  = toPtr(oldnum & ~MARKBITS);
        oldmark = oldnum & MARKBITS;
        return res;
      }

      /// String compare_and_set based on old and new pairs of pointer and mark values.
      /// \note Comparable to AtomicMarkableReference.CompareAndSet.
      bool compare_exchange_strong( T*&               oldptr,
                                    T*                newptr,
                                    mark_type&        oldmark,
                                    mark_type         newmark,
                                    std::memory_order succ = std::memory_order_seq_cst,
                                    std::memory_order fail = std::memory_order_seq_cst
                                  )
      {
        assert(  (oldmark < MAX_MARKED)
              && ((toNum(oldptr) & MARKBITS) == 0)
              && (newmark < MAX_MARKED)
              && ((toNum(newptr) & MARKBITS) == 0)
              );

        numptr_type       oldnum = reinterpret_cast<numptr_type>(oldptr) | oldmark;
        const numptr_type newnum = reinterpret_cast<numptr_type>(newptr) | newmark;
        const bool     res    =  val.compare_exchange_strong(oldnum, newnum, succ, fail);

        oldptr  = toPtr(oldnum & ~MARKBITS);
        oldmark = oldnum & MARKBITS;
        return res;
      }

      /// Reads a pointer's mark state.
      /// \note Comparable to AtomicMarkableReference.isMarked
      mark_type mark(std::memory_order order = std::memory_order_seq_cst) const
      {
        return static_cast<mark_type>(val.load(order) & MARKBITS);
      }

      /// Returns the current state of the pointer as pointer and mark pair.
      /// comparable to AtomicMarkableReference.get
      state_type
      state(std::memory_order order = std::memory_order_seq_cst) const
      {
        numptr_type v = val.load(order);

        return std::make_pair(toPtr(v & ~MARKBITS), static_cast<mark_type>(v & MARKBITS));
      }

      /// sequentially consistent pointer dereference
      T& operator*()
      {
        return *load();
      }

      /// sequentially consistent pointer load
      T* operator->()
      {
        return load();
      }

    private:
      std::atomic<numptr_type> val;
  };
} // end namespace uab
#endif /* _ATOMICUTIL_HPP */
