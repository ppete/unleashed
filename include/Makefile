
# bitutil.hpp circularqueue.hpp  vector.hpp
SOURCES=pmemory.hpp skiplist.hpp atomicutil.hpp

.PHONY: dox clean
.PRECIOUS: Doxyfile

default: dox

dox: Doxyfile.in $(SOURCES)
	mkdir -p ./doc
	sed 's/UABDOXDIR/.\/doc/g' $< >Doxyfile
	sed -i 's/UABSRC/$(SOURCES)/g' Doxyfile
	doxygen

clean:
	rm -rf ./doc

