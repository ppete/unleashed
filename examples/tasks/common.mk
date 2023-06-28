
##
## expects SOURCES variable be set to list of files to be compiled

# extract name of the main source file
CODE:=$(basename $(notdir $(firstword $(SOURCES))))

##
## compiler name

COMP:=$(notdir $(CXX))

##
## output directory

ifeq ($(COMPDIR),)
  THEHOST:=$(shell hostname -s | sed -e s/[0-9]*//g)
  COMPDIR:=$(UCL_HOME)/examples/tasks/tmp/$(THEHOST)
  dummy := $(shell mkdir -p $(COMPDIR))
  $(info *** COMPDIR set to $(COMPDIR))
endif


##
## Makefile TARGETS

.PHONY: default clean

TARGETS = $(COMPDIR)/$(CODE)-ucl-$(COMP).bin 
          

ifeq ($(TBB_ENABLED),1)
  TARGETS += $(COMPDIR)/$(CODE)-tbb-$(COMP).bin
endif

ifeq ($(CILK_ENABLED),1)
  TARGETS += $(COMPDIR)/$(CODE)-cilk-$(COMP).bin
endif

ifeq ($(QTHREADS_ENABLED),1)
  TARGETS += $(COMPDIR)/$(CODE)-qthreads-$(COMP).bin
endif

ifeq ($(OPENMP_ENABLED),1)
  TARGETS += $(COMPDIR)/$(CODE)-omp-$(COMP).bin

  ifeq ($(WOMP_ENABLED),1)
    TARGETS += $(COMPDIR)/$(CODE)-womp-$(COMP).bin
  endif
endif

ifeq ($(SEQ_ENABLED),1)
  TARGETS += $(COMPDIR)/$(CODE)-seq-$(COMP).bin

  ifeq ($(WOMP_ENABLED),1)
    OPENMPV:=-DWOMP_VERSION=1
  else
    OPENMPV:=-DOMP_VERSION=1
  endif
endif


CXXVERSION ?= -std=c++11

default: $(TARGETS)

CXXFLAGS = $(CXXVERSION) $(WARNFLAG) $(OPTFLAG) $(CPUARCH) $(DBGFLAG)

$(COMPDIR)/$(CODE)-seq-$(COMP).bin: $(SOURCES)
	$(CXX) $(CXXFLAGS) $(OPENMPV) -I$(UCL_HOME)/include $(SOURCES) $(LINKATOMIC) -o $@

$(COMPDIR)/$(CODE)-ucl-$(COMP).bin: $(SOURCES) $(UCL_HOME)/include/ucl/task.hpp $(UCL_HOME)/include/ucl/task-pool-x.hpp
	$(CXX) $(CXXFLAGS) $(THREADFLAG) $(UCLFLAG) -DUCL_VERSION=1 -I$(UCL_HOME)/include $(SOURCES) $(LINKATOMIC) -o $@

$(COMPDIR)/$(CODE)-omp-$(COMP).bin: $(SOURCES)
	$(CXX) $(CXXFLAGS) $(OMPFLAG) -DOMP_VERSION=1 -I$(UCL_HOME)/include $(SOURCES) $(LINKATOMIC) -o $@

$(COMPDIR)/$(CODE)-womp-$(COMP).bin: $(SOURCES)
	$(CXX) $(CXXFLAGS) $(OMPFLAG) -DWOMP_VERSION=1 -I$(UCL_HOME)/include $(SOURCES) $(LINKATOMIC) -o $@

$(COMPDIR)/$(CODE)-cilk-$(COMP).bin: $(SOURCES)
	$(CXX) $(CXXFLAGS) $(CILKFLAG) -DCILK_VERSION=1 -I$(UCL_HOME)/include $(SOURCES) $(LINKATOMIC) -o $@

$(COMPDIR)/$(CODE)-tbb-$(COMP).bin: $(SOURCES)
	$(CXX) $(CXXFLAGS) $(THREADFLAG) -DTBB_VERSION=1 -I$(UCL_HOME)/include -I$(TBB_HOME)/include $(SOURCES) $(LINKATOMIC) -L$(TBB_HOME)/lib -ltbb -ltbbmalloc -o $@

$(COMPDIR)/$(CODE)-qthreads-$(COMP).bin: $(SOURCES)
	$(CXX) $(CXXFLAGS) $(THREADFLAG) -DQTHREADS_VERSION=1 -I$(UCL_HOME)/include -I$(QTHREADS_HOME)/include $(SOURCES) $(LINKATOMIC) -L$(QTHREADS_HOME)/lib -lqthread -o $@


clean:
	rm -f $(COMPDIR)/$(CODE)*.bin
