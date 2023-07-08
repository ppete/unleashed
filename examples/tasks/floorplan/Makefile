include ../../make-util/discover-toolset.mk
include ../../make-util/discover-task-models.mk

ifeq ($(CXX17_ENABLEDX),1)	
	CXXVERSION:=-std=c++17
	
	# derived from BOTS (include its original version)
	WOMP_ENABLED:=$(OPENMP_ENABLED)
	
	# source files
	SOURCES:=floorplan.cc
	
	include ../common.mk
else
  default:
  $(info Floorplan benchmark not compiled- requires C++17 support)	
endif	
	
