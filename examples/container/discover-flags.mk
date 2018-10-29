

ifeq ($(BLAZE_HOME),)
  export BLAZE_HOME=../..
  $(info *** BLAZE_HOME set to default value)
endif


##
## auto set instruction set to native
## with Clang and Power architecture, this does not work (use INSTRSET instead)

ifeq ($(INSTRSET),)
	SUCCESS=$(shell $(CXX) -std=c++11 -O0 -march=native $(BLAZE_HOME)/examples/common/test-hello.cc -o $(OUTPUTDIR)/test-hello.bin; echo $$?)

	ifeq ($(SUCCESS),0)
    export INSTRSET= -march=native
    $(info *** using $(INSTRSET))
  endif
endif

ifeq ($(INSTRSET),)
	SUCCESS=$(shell $(CXX) -std=c++11 -O0 -mcpu=native $(BLAZE_HOME)/examples/common/test-hello.cc -o $(OUTPUTDIR)/test-hello.bin; echo $$?)

	ifeq ($(SUCCESS),0)
    export INSTRSET= -mcpu=native
    $(info *** using $(THREADFLAG))
  endif
endif


##
## auto set threading approach (sun c++ vs everyone else)
ifeq ($(THREADFLAG),)
  SUCCESS=$(shell $(CXX) -std=c++11 -O0 -pthread $(BLAZE_HOME)/examples/common/test-hello.cc -o $(OUTPUTDIR)/test-hello.bin; echo $$?)

  ifeq ($(SUCCESS),0)
    export THREADFLAG= -pthread
    $(info *** using threads $(INSTRSET))
  endif
endif

ifeq ($(THREADFLAG),)
  SUCCESS=$(shell $(CXX) -std=c++11 -O0 -mt $(BLAZE_HOME)/examples/common/test-hello.cc -o $(OUTPUTDIR)/test-hello.bin; echo $$?)

  ifeq ($(SUCCESS),0)
    export THREADFLAG= -mt
    $(info *** using threads $(THREADFLAG))
  endif
endif


##
## auto set HTM, if TEST_HTM is enabled

ifeq ($(TEST_HTM),)
  TXFLAGS=
else
  DEFINES+= -DHTM_ENABLED=1

  ifeq ($(THISARCH),)
    $(error THISARCH not set (needs to be INTEL or POWER8))
  endif

  # try INTEL flag
  ifeq ($(TXFLAGS),)
    SUCCESS=$(shell $(CXX) -std=c++11 -O0 -mrtm $(BLAZE_HOME)/examples/common/test-hello.cc -o $(OUTPUTDIR)/test-hello.bin; echo $$?)

    ifeq ($(SUCCESS),0)
      export TXFLAGS= -mrtm
      $(info *** using threads $(TXFLAGS))
    endif
  endif

  # try IBM POWER flag
  ifeq ($(TXFLAGS),)
    SUCCESS=$(shell $(CXX) -std=c++11 -O0 -mhtm $(BLAZE_HOME)/examples/common/test-hello.cc -o $(OUTPUTDIR)/test-hello.bin; echo $$?)

    ifeq ($(SUCCESS),0)
      export TXFLAGS= -mhtm
      $(info *** using threads $(TXFLAGS))
    endif
  endif

  ifeq ($(TXFLAGS),)
    $(error *** using to identify hardware transactional memory. try setting TXFLAGS)
  endif
endif
