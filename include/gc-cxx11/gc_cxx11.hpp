#ifndef CXX11_GC_ALLOCATOR_HPP
#define CXX11_GC_ALLOCATOR_HPP

#include "gc_allocator.h"

#ifdef BLAZE_GC_CXX11_THREAD_CONTEXT
#error "GC_CXX_THREAD_CONTEXT already defined. Check include order."
#endif /* BLAZE_GC_CXX11_THREAD_CONTEXT */

#define BLAZE_GC_CXX11_THREAD_CONTEXT 1

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

  private:
    void* operator new(size_t size) = delete;
};

#endif /* CXX11_GC_ALLOCATOR_HPP */
