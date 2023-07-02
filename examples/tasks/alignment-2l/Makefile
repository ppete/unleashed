
include ../../make-util/discover-toolset.mk
include ../../make-util/discover-task-models.mk

# derived from BOTS (include its original version)
WOMP_ENABLED:=1

# yield segmentation faults on stack allocations with 5000 elements
#   observed on Intel dual socket, gcc 7.3
QTHREADS_ENABLED:=0

# source files
SOURCES:=alignment.cc

include ../common.mk
