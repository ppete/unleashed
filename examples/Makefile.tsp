ifeq ($(OUTPUTDIR),)
	OUTPUTDIR=./tmp
endif

ifeq ($(BLAZE_HOME),)
	BLAZE_HOME=..
endif

ifneq ($(NUMTHREADS),)
	EXTRAFLAGS=-DNUMTHREADS=$(NUMTHREADS)
	SUFFIX=.$(NUMTHREADS)
else
	EXTRAFLAGS=
	SUFFIX=
endif

ifneq ($(NUM_NODES),)
  EXTRAFLAGS+=-DNUM_NODES=$(NUM_NODES)
endif

ifeq ($(TBB_HOME),)
  TBB_HOME=~/../../../usr
endif

#GLOB=
#GLOB=-DZEROGLOBALS=1 # turns of global result reduction for TBB
GLOB= #ADDED

.PHONY: default clean

default: \
	 $(OUTPUTDIR)/tsp-omp$(SUFFIX).bin \
	 $(OUTPUTDIR)/tsp-blaze$(SUFFIX).bin \
	 $(OUTPUTDIR)/tsp-tbb$(SUFFIX).bin


#~ CXXFLAGS = -ggdb -std=c++11 -Wall -Wextra -pedantic -O2 -DNDEBUG=1 -march=native

CXXFLAGS = -std=c++11 -Wall -Wextra -pedantic -O2 -DNDEBUG=1 -march=native $(EXTRAFLAGS) #ADDED -fopenmp

$(OUTPUTDIR)/tsp-blaze$(SUFFIX).bin: tsp.cc Makefile.tsp
	$(CXX) $(CXXFLAGS) -pthread -DBLAZE_VERSION=1 -I$(BLAZE_HOME)/include $< -o $@

$(OUTPUTDIR)/tsp-omp$(SUFFIX).bin: tsp.cc Makefile.tsp
	$(CXX) $(CXXFLAGS) -fopenmp -DOMP_VERSION=1 -I$(BLAZE_HOME)/include $< -o $@

$(OUTPUTDIR)/tsp-tbb$(SUFFIX).bin: tsp.cc Makefile.tsp
	$(CXX) $(CXXFLAGS) -pthread -DTBB_VERSION=1 -I$(BLAZE_HOME)/include -I$(TBB_HOME)/include -L$(TBB_HOME)/lib $(GLOB) $< -ltbb -ltbbmalloc -o $@

clean:
	rm -rf $(OUTPUTDIR)/*
