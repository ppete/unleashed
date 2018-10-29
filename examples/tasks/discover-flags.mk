##
## set directories if not set

ifeq ($(TBB_HOME),)
  export TBB_HOME=/usr
  $(info *** TBB_HOME set to default value)
endif

ifeq ($(BLAZE_HOME),)
  export BLAZE_HOME=../../..
  $(info *** BLAZE_HOME set to default value)
endif

ifeq ($(QTHREADS_HOME),)
  export QTHREADS_HOME=/usr
  $(info *** QTHREADS_HOME set to default value)
endif


ifeq ($(OUTPUTDIR),)
  export OUTPUTDIR=$(BLAZE_HOME)/examples/tasks/tmp
  $(info *** OUTPUTDIR set to $(OUTPUTDIR))
endif


##
## set suffix

ifneq ($(NUMTHREADS),)
  export EXTRAFLAGS=-DNUMTHREADS=$(NUMTHREADS)
  export SUFFIX=-$(NUMTHREADS)
	ifeq ($(TEST_HTM),)
	  export BLZSUFFIX=-$(NUMTHREADS)
	else
	  export BLZSUFFIX=-$(NUMTHREADS)-tx
	endif
else
  export EXTRAFLAGS=
  export SUFFIX=
endif


##
## set compiler flags

export HTMFLAGS=


##
## use Blaze/HTM version

ifeq ($(TEST_HTM),1)
  ifeq ($(THISARCH),)
    $(error THISARCH not set (needs to be INTEL or POWER))
  endif

	ifeq ($(THISARCH), INTEL)
		HTMFLAGS+= -mrtm
	endif

	ifeq ($(THISARCH), POWER)
		HTMFLAGS+= -mhtm
	endif

	HTMFLAGS+= -DHTM_ENABLED=1
endif

##
## set problem size

ifneq ($(PROBLEM_SIZE),)
	EXTRAFLAGS+= -DPROBLEM_SIZE=$(PROBLEM_SIZE)
endif


##
## auto enabling programming models

ifeq ($(CILK_ENABLED),)
  SUCCESS=$(shell $(CXX) -std=c++11 -O0 -fcilkplus -lcilkrts $(BLAZE_HOME)/examples/tasks/common/test-cilk.cc -o $(BLAZE_HOME)/examples/tasks/common/test-cilk.bin; echo $$?)

  ifeq ($(SUCCESS),0)
    export CILK_ENABLED=1
  else
    export CILK_ENABLED=0
    $(info *** CILK disabled)
  endif
endif

ifeq ($(QTHREADS_ENABLED),)
  SUCCESS=$(shell $(CXX) -std=c++11 -O0 -I$(QTHREADS_HOME)/include $(BLAZE_HOME)/examples/tasks/common/test-qthreads.cc -L$(QTHREADS_HOME)/lib -lqthread -o $(BLAZE_HOME)/examples/tasks/common/test-qthreads.bin; echo $$?)

  ifeq ($(SUCCESS),0)
    export QTHREADS_ENABLED=1
  else
    export QTHREADS_ENABLED=0
    $(info *** QTHREADS disabled)
  endif
endif


ifeq ($(TBB_ENABLED),)
  SUCCESS=$(shell $(CXX) -std=c++11 -O0 -pthread -I$(TBB_HOME)/include $(BLAZE_HOME)/examples/tasks/common/test-tbb.cc -latomic -L$(TBB_HOME)/lib -ltbb -ltbbmalloc -o $(BLAZE_HOME)/examples/tasks/common/test-tbb.bin; echo $$?)

  ifeq ($(SUCCESS),0)
    export TBB_ENABLED=1
  else
    export TBB_ENABLED=0
    $(info *** TBB disabled)
  endif
endif

##
## auto set instruction set to native

ifeq ($(INSTRSET),)
	SUCCESS=$(shell $(CXX) -std=c++11 -O0 -march=native $(BLAZE_HOME)/examples/tasks/common/test-hello.cc -o $(BLAZE_HOME)/examples/tasks/common/test-hello.bin; echo $$?)

	ifeq ($(SUCCESS),0)
    export INSTRSET= -march=native
    $(info *** using $(INSTRSET))
  endif
endif

ifeq ($(INSTRSET),)
	SUCCESS=$(shell $(CXX) -std=c++11 -O0 -mcpu=native $(BLAZE_HOME)/examples/tasks/common/test-hello.cc -o $(BLAZE_HOME)/examples/tasks/common/test-hello.bin; echo $$?)

	ifeq ($(SUCCESS),0)
    export INSTRSET= -mcpu=native
    $(info *** using $(THREADFLAG))
  endif
endif


##
## auto set threading approach (sun c++ vs everyone else)
ifeq ($(THREADFLAG),)
  SUCCESS=$(shell $(CXX) -std=c++11 -O0 -pthread $(BLAZE_HOME)/examples/tasks/common/test-hello.cc -o $(BLAZE_HOME)/examples/tasks/common/test-hello.bin; echo $$?)

  ifeq ($(SUCCESS),0)
    export THREADFLAG= -pthread
    $(info *** using threads $(INSTRSET))
  endif
endif

ifeq ($(THREADFLAG),)
  SUCCESS=$(shell $(CXX) -std=c++11 -O0 -mt $(BLAZE_HOME)/examples/tasks/common/test-hello.cc -o $(BLAZE_HOME)/examples/tasks/common/test-hello.bin; echo $$?)

  ifeq ($(SUCCESS),0)
    export THREADFLAG= -mt
    $(info *** using threads $(THREADFLAG))
  endif
endif
