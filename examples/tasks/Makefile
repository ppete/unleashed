
ifeq ($(NUMTHREADS),)
  NUMTHREADS=20
endif

include ../make-util/discover-toolset.mk
include ../make-util/discover-task-models.mk

.PHONY: default clean

default:
	mkdir -p ./tmp
	cd alignment-2l; $(MAKE) OUTPUTDIR=../tmp NUMTHREADS=$(NUMTHREADS)
	cd dot-product; $(MAKE) OUTPUTDIR=../tmp NUMTHREADS=$(NUMTHREADS)
	cd fib; $(MAKE) OUTPUTDIR=../tmp NUMTHREADS=$(NUMTHREADS)
	cd floorplan; $(MAKE) OUTPUTDIR=../tmp NUMTHREADS=$(NUMTHREADS)
	cd knapsack; $(MAKE) OUTPUTDIR=../tmp NUMTHREADS=$(NUMTHREADS)
	cd nqueens-0; $(MAKE) OUTPUTDIR=../tmp NUMTHREADS=$(NUMTHREADS)
	cd pi; $(MAKE) OUTPUTDIR=../tmp NUMTHREADS=$(NUMTHREADS)
	cd tsp; $(MAKE) OUTPUTDIR=../tmp NUMTHREADS=$(NUMTHREADS)
	cd uts; $(MAKE) OUTPUTDIR=../tmp NUMTHREADS=$(NUMTHREADS)

clean:
	cd alignment-2l; $(MAKE) OUTPUTDIR=../tmp clean
	cd dot-product; $(MAKE) OUTPUTDIR=../tmp clean
	cd fib; $(MAKE) OUTPUTDIR=../tmp clean
	cd floorplan; $(MAKE) OUTPUTDIR=../tmp clean
	cd knapsack; $(MAKE) OUTPUTDIR=../tmp clean
	cd nqueens-0; $(MAKE) OUTPUTDIR=../tmp clean
	cd pi; $(MAKE) OUTPUTDIR=../tmp clean
	cd tsp; $(MAKE) OUTPUTDIR=../tmp clean
	cd uts; $(MAKE) OUTPUTDIR=../tmp clean
