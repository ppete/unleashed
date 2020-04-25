
include ../make-util/discover-toolset.mk
include ../make-util/discover-task-models.mk

.PHONY: default clean

default:
	mkdir -p ./tmp
	cd pi; $(MAKE)
	cd alignment-2l; $(MAKE)
	cd dot-product; $(MAKE)
	cd fib; $(MAKE)
	cd floorplan; $(MAKE)
	cd knapsack; $(MAKE)
	cd nqueens-0; $(MAKE)
	cd qsort; $(MAKE)
	cd tsp; $(MAKE)
	cd uts; $(MAKE)

clean:
	cd alignment-2l; $(MAKE)  clean
	cd dot-product; $(MAKE)  clean
	cd fib; $(MAKE)  clean
	cd floorplan; $(MAKE)  clean
	cd knapsack; $(MAKE)  clean
	cd msort; $(MAKE)  clean
	cd nqueens-0; $(MAKE)  clean
	cd pi; $(MAKE)  clean
	cd tsp; $(MAKE)  clean
	cd uts; $(MAKE)  clean
