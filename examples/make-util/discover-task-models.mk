##
## set directories if not set

#~ include ../make-util/discover-toolset.mk

ifeq ($(TBB_HOME),)
  export TBB_HOME=/usr
  $(info *** TBB_HOME set to default value)
endif

ifeq ($(QTHREADS_HOME),)
  export QTHREADS_HOME=/usr
  $(info *** QTHREADS_HOME set to default value)
endif

##
## auto enabling programming models

ifeq ($(CILK_ENABLED),)
  SUCCESS=$(shell $(CXX) -std=c++11 -O0 -fcilkplus -lcilkrts $(UCL_HOME)/examples/make-util/test-cilk.cc -o /dev/null; echo $$?)

  ifeq ($(SUCCESS),0)
    export CILK_ENABLED=1
  else
    export CILK_ENABLED=0
    $(info *** CILK disabled)
  endif
endif

ifeq ($(QTHREADS_ENABLED),)
  SUCCESS=$(shell $(CXX) -std=c++11 -O0 -I$(QTHREADS_HOME)/include $(UCL_HOME)/examples/make-util/test-qthreads.cc -L$(QTHREADS_HOME)/lib -lqthread -o /dev/null; echo $$?)

  ifeq ($(SUCCESS),0)
    export QTHREADS_ENABLED=1
  else
    export QTHREADS_ENABLED=0
    $(info *** QTHREADS disabled)
  endif
endif


ifeq ($(TBB_ENABLED),)
  SUCCESS=$(shell $(CXX) -std=c++11 -O0 -I$(TBB_HOME)/include $(UCL_HOME)/examples/make-util/test-tbb.cc -latomic -L$(TBB_HOME)/lib -ltbb -ltbbmalloc -o /dev/null; echo $$?)
#~   $(info "$(CXX) -std=c++11 -O0 -I$(TBB_HOME)/include $(UCL_HOME)/examples/make-util/test-tbb.cc -latomic -L$(TBB_HOME)/lib -ltbb -ltbbmalloc -o $(UCL_HOME)/examples/make-util/test-tbb.bin")

  ifeq ($(SUCCESS),0)
    export TBB_ENABLED=1
  else
    export TBB_ENABLED=0
    $(info *** TBB disabled)
  endif
endif
