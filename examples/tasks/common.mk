
##
## expects SOURCES variable be set to list of files to be compiled

# extract name of the main source file
CODE=$(basename $(notdir $(firstword $(SOURCES))))

##
## compiler name

COMP=$(notdir $(CXX))

##
## output directory

ifeq ($(COMPDIR),)
  export COMPDIR=$(UCL_HOME)/examples/tasks/tmp
  $(info *** COMPDIR set to $(COMPDIR))
endif


##
## Makefile TARGETS

.PHONY: default clean

TARGETS = $(COMPDIR)/$(CODE)-ucl-$(COMP).bin \
          $(COMPDIR)/$(CODE)-omp-$(COMP).bin

ifeq ($(TBB_ENABLED),1)
  TARGETS += $(COMPDIR)/$(CODE)-tbb-$(COMP).bin
endif

ifeq ($(CILK_ENABLED),1)
  TARGETS += $(COMPDIR)/$(CODE)-cilk-$(COMP).bin
endif

ifeq ($(QTHREADS_ENABLED),1)
  TARGETS += $(COMPDIR)/$(CODE)-qthreads-$(COMP).bin
endif

ifeq ($(WOMP_ENABLED),1)
  TARGETS += $(COMPDIR)/$(CODE)-womp-$(COMP).bin
endif

ifeq ($(SEQ_ENABLED),1)
  TARGETS += $(COMPDIR)/$(CODE)-seq-$(COMP).bin

  ifeq ($(WOMP_ENABLED),1)
    OPENMPV=-DWOMP_VERSION=1
  else
    OPENMPV=-DOMP_VERSION=1
  endif
endif

ifeq ($(DBGFLAG),)
  DBGFLAG:=-DNDEBUG=1
else
  OPTFLAG:=-O0
  $(info *** changed OPTFLAG to -Og)
endif

default: $(TARGETS)

CXXFLAGS = -std=c++11 $(WARNFLAG) $(OPTFLAG) $(CPUARCH) $(DBGFLAG)

$(COMPDIR)/$(CODE)-seq-$(COMP).bin: $(SOURCES)
	$(CXX) $(CXXFLAGS) $(OPENMPV) -I$(UCL_HOME)/include $(SOURCES) $(LINKATOMIC) -o $@

$(COMPDIR)/$(CODE)-ucl-$(COMP).bin: $(SOURCES) $(UCL_HOME)/include/ucl/task.hpp $(UCL_HOME)/include/ucl/task-pool-x.hpp
	$(CXX) $(CXXFLAGS) $(THREADFLAG) -DUCL_VERSION=1 -I$(UCL_HOME)/include $(SOURCES) $(LINKATOMIC) -o $@

$(COMPDIR)/$(CODE)-htmucl-$(COMP).bin: $(SOURCES) $(UCL_HOME)/include/ucl/task.hpp
	$(CXX) $(CXXFLAGS) $(THREADFLAG) $(HTMFLAG) -DUCL_VERSION=1 -I$(UCL_HOME)/include $(SOURCES) $(LINKATOMIC) -o $@

$(COMPDIR)/$(CODE)-omp-$(COMP).bin: $(SOURCES)
	$(CXX) $(CXXFLAGS) $(OMPFLAG) -DOMP_VERSION=1 -I$(UCL_HOME)/include $(SOURCES) $(LINKATOMIC) -o $@

$(COMPDIR)/$(CODE)-womp-$(COMP).bin: $(SOURCES)
	$(CXX) $(CXXFLAGS) $(OMPFLAG) -DWOMP_VERSION=1 -I$(UCL_HOME)/include $(SOURCES) $(LINKATOMIC) -o $@

$(COMPDIR)/$(CODE)-cilk-$(COMP).bin: $(SOURCES)
	$(CXX) $(CXXFLAGS) $(CILKFLAG) -DCILK_VERSION=1 -I$(UCL_HOME)/include $(SOURCES) $(LINKATOMIC) -o $@

$(COMPDIR)/$(CODE)-tbb-$(COMP).bin: $(SOURCES)
	$(CXX) $(CXXFLAGS) $(THREADFLAG) -DTBB_VERSION=1 -I$(UCL_HOME)/include -I$(TBB_HOME)/include $(SOURCES) -latomic -L$(TBB_HOME)/lib -ltbb -ltbbmalloc -o $@

$(COMPDIR)/$(CODE)-qthreads-$(COMP).bin: $(SOURCES)
	$(CXX) $(CXXFLAGS) $(THREADFLAG) -DQTHREADS_VERSION=1 -I$(UCL_HOME)/include -I$(QTHREADS_HOME)/include $(SOURCES) $(LINKATOMIC) -L$(QTHREADS_HOME)/lib -lqthread -o $@


clean:
	rm -f $(COMPDIR)/$(CODE)*.bin
