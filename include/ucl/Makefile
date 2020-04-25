
# bitutil.hpp circularqueue.hpp  vector.hpp
SOURCES=pmemory.hpp skiplist.hpp atomicutil.hpp task.hpp queue.hpp stack.hpp spinlock.hpp \
        htm.hpp htm-skiplist.hpp htm-memory.hpp bitutil.hpp archmodel.hpp

.PHONY: dox clean
.PRECIOUS: Doxyfile

default: dox

dox: Doxyfile.in $(SOURCES)
	mkdir -p /tmp/doc
	sed 's/UCLDOXDIR/.\/doc/g' $< >Doxyfile
	sed -i 's/UCLSRC/$(SOURCES)/g' Doxyfile
	doxygen
	mv doc/html/index.html doc/html/main.html
	cp doc/html/files.html doc/html/index.html

clean:
	rm -rf ./doc
