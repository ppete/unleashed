
##
## expects SOURCES variable be set to list of files to be compiled

# extract name of the main source file
CODE=$(basename $(notdir $(firstword $(SOURCES))))

##
## set directories if not set

ifeq ($(TBB_HOME),)
  TBB_HOME=/usr
  $(info *** TBB_HOME set to default value)
endif

ifeq ($(BLAZE_HOME),)
  BLAZE_HOME=../../..
  $(info *** BLAZE_HOME set to default value)
endif

ifeq ($(QTHREADS_HOME),)
  QTHREADS_HOME=/usr
  $(info *** QTHREADS_HOME set to default value)
endif


ifeq ($(OUTPUTDIR),)
  OUTPUTDIR=$(BLAZE_HOME)/examples/tasks/tmp
  $(info *** OUTPUTDIR set to $(OUTPUTDIR))
endif


##
## set suffix

ifneq ($(NUMTHREADS),)
  EXTRAFLAGS=-DNUMTHREADS=$(NUMTHREADS)
  SUFFIX=-$(NUMTHREADS)
	ifeq ($(TEST_HTM),)
	  BLZSUFFIX=-$(NUMTHREADS)
	else
	  BLZSUFFIX=-$(NUMTHREADS)-tx
	endif
else
  EXTRAFLAGS=
  SUFFIX=
endif


##
## set compiler flags

OMPLINKFLAGS=
HTMFLAGS=


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
    CILK_ENABLED=1
  else
    $(info *** CILK disabled)
  endif
endif

ifeq ($(QTHREADS_ENABLED),)
  SUCCESS=$(shell $(CXX) -std=c++11 -O0 -I$(QTHREADS_HOME)/include $(BLAZE_HOME)/examples/tasks/common/test-qthreads.cc -L$(QTHREADS_HOME)/lib -lqthread -o $(BLAZE_HOME)/examples/tasks/common/test-qthreads.bin; echo $$?)

  ifeq ($(SUCCESS),0)
    QTHREADS_ENABLED=1
  else
    $(info *** QTHREADS disabled)
  endif
endif


ifeq ($(TBB_ENABLED),)
  SUCCESS=$(shell $(CXX) -std=c++11 -O0 -pthread -I$(TBB_HOME)/include $(BLAZE_HOME)/examples/tasks/common/test-tbb.cc -latomic -L$(TBB_HOME)/lib -ltbb -ltbbmalloc -o $(BLAZE_HOME)/examples/tasks/common/test-tbb.bin; echo $$?)

  ifeq ($(SUCCESS),0)
    TBB_ENABLED=1
  else
    $(info *** TBB disabled)
  endif
endif

##
## auto set instruction set to native

ifeq ($(INSTRSET),)
	SUCCESS=$(shell $(CXX) -std=c++11 -O0 -march=native $(BLAZE_HOME)/examples/tasks/common/test-hello.cc -o $(BLAZE_HOME)/examples/tasks/common/test-hello.bin; echo $$?)

	ifeq ($(SUCCESS),0)
    INSTRSET= -march=native
    $(info *** using $(INSTRSET))
  endif
endif

ifeq ($(INSTRSET),)
	SUCCESS=$(shell $(CXX) -std=c++11 -O0 -mcpu=native $(BLAZE_HOME)/examples/tasks/common/test-hello.cc -o $(BLAZE_HOME)/examples/tasks/common/test-hello.bin; echo $$?)

	ifeq ($(SUCCESS),0)
    INSTRSET= -mcpu=native
    $(info *** using $(INSTRSET))
  endif
endif


##
## compiler name

COMP=$(notdir $(CXX))


##
## Makefile TARGETS

.PHONY: default clean

TARGETS = $(OUTPUTDIR)/$(CODE)-blaze-$(COMP)$(BLZSUFFIX).bin \
          $(OUTPUTDIR)/$(CODE)-omp-$(COMP)$(SUFFIX).bin

ifeq ($(TBB_ENABLED),1)
  TARGETS += $(OUTPUTDIR)/$(CODE)-tbb-$(COMP)$(SUFFIX).bin
endif

ifeq ($(CILK_ENABLED),1)
  TARGETS += $(OUTPUTDIR)/$(CODE)-cilk-$(COMP)$(SUFFIX).bin
endif

ifeq ($(QTHREADS_ENABLED),1)
  TARGETS += $(OUTPUTDIR)/$(CODE)-qthreads-$(COMP)$(SUFFIX).bin
endif

ifeq ($(BOTS_ORIGINAL),1)
  TARGETS += $(OUTPUTDIR)/$(CODE)-bots-$(COMP)$(SUFFIX).bin
endif


default: $(TARGETS)

#~ CXXFLAGS = -std=c++11 -ggdb -pg -Wall -Wextra -pedantic -O2 $(INSTRSET) $(EXTRAFLAGS)

#~ CXXFLAGS = -std=c++11 -ggdb -Wall -Wextra -pedantic -Og $(INSTRSET) $(EXTRAFLAGS)

CXXFLAGS = -std=c++11 -Wall -Wextra -pedantic -O2 $(INSTRSET) -DNDEBUG=1 $(EXTRAFLAGS)

$(OUTPUTDIR)/$(CODE)-blaze-$(COMP)$(BLZSUFFIX).bin: $(SOURCES)
	$(CXX) $(CXXFLAGS) $(HTMFLAGS) -pthread -DBLAZE_VERSION=1 -I$(BLAZE_HOME)/include $(SOURCES) $(LINKATOMIC) -o $@

$(OUTPUTDIR)/$(CODE)-omp-$(COMP)$(SUFFIX).bin: $(SOURCES)
	$(CXX) $(CXXFLAGS) -fopenmp -DOMP_VERSION=1 -I$(BLAZE_HOME)/include $(SOURCES) $(OMPLINKFLAGS) $(LINKATOMIC) -o $@

$(OUTPUTDIR)/$(CODE)-bots-$(COMP)$(SUFFIX).bin: $(SOURCES)
	$(CXX) $(CXXFLAGS) -fopenmp -DBOTS_VERSION=1 -I$(BLAZE_HOME)/include $(SOURCES) $(OMPLINKFLAGS) $(LINKATOMIC) -o $@

$(OUTPUTDIR)/$(CODE)-cilk-$(COMP)$(SUFFIX).bin: $(SOURCES)
	$(CXX) $(CXXFLAGS) -fcilkplus -lcilkrts -DCILK_VERSION=1 -I$(BLAZE_HOME)/include $(SOURCES) $(LINKATOMIC) -o $@

$(OUTPUTDIR)/$(CODE)-tbb-$(COMP)$(SUFFIX).bin: $(SOURCES)
	$(CXX) $(CXXFLAGS) -pthread -DTBB_VERSION=1 -I$(BLAZE_HOME)/include -I$(TBB_HOME)/include $(SOURCES) -latomic -L$(TBB_HOME)/lib -ltbb -ltbbmalloc -o $@

$(OUTPUTDIR)/$(CODE)-qthreads-$(COMP)$(SUFFIX).bin: $(SOURCES)
	$(CXX) $(CXXFLAGS) -pthread -DQTHREADS_VERSION=1 -I$(BLAZE_HOME)/include -I$(QTHREADS_HOME)/include $(SOURCES) $(LINKATOMIC) -L$(QTHREADS_HOME)/lib -lqthread -o $@


clean:
	rm -f $(OUTPUTDIR)/$(CODE)*bin
	rm -f $(BLAZE_HOME)/examples/tasks/common/test*.bin
