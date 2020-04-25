
##
## discover UCL HOME

ifeq ($(UCL_HOME),)
  ALOCATION:=$(lastword $(MAKEFILE_LIST))
  BLOCATION:=$(realpath $(ALOCATION))
  MYLOCATION:=$(subst discover-toolset.mk,,$(BLOCATION))
#~   MYLOCATION=$(subst common.mk,,$(CLOCATION))

  export UCL_HOME:=$(MYLOCATION)/../..
endif


##
## auto set system

ifeq ($(TARGETARCH),)
  THISSYSTEM:=$(shell uname -a)

  ifneq (,$(findstring power,$(THISSYSTEM)))
    TARGETARCH:=POWER
  endif

  ifneq (,$(findstring ppc,$(THISSYSTEM)))
    TARGETARCH:=POWER
  endif

  ifneq (,$(findstring 386,$(THISSYSTEM)))
    TARGETARCH:=i386
  endif

  ifneq (,$(findstring x86_64,$(THISSYSTEM)))
    TARGETARCH:=x86_64
  endif

  ifneq (,$(findstring amd,$(THISSYSTEM)))
    TARGETARCH:=AMD
  endif

  ifneq (,$(findstring sparc,$(THISSYSTEM)))
    TARGETARCH:=SPARC
  endif
endif


##
## auto set toolset and flags
ifeq ($(TOOLSET),)
  ifeq ($(COMPILERVERSION),)
    SUCCESS:=$(shell $(CXX) -V -c >/dev/null 2>&1; echo $$?)

    ifeq ($(SUCCESS),0)
      COMPILERVERSION:=$(shell $(CXX) -V 2>&1)
    endif
  endif

  ifeq ($(COMPILERVERSION),)
    SUCCESS:=$(shell $(CXX) --version -c >/dev/null 2>&1; echo $$?)

    ifeq ($(SUCCESS),0)
      COMPILERVERSION:=$(shell $(CXX) --version)
    endif
  endif
endif


ifeq ($(TOOLSET),)
  ifneq (,$(findstring Free Software Foundation,$(COMPILERVERSION)))
    export TOOLSET    := GCC
    export OPTFLAG    ?= -O2
    export DBGOPTFLAG ?= -O0
    export THREADFLAG ?= -pthread
    export WARNFLAG   ?= -Wall -Wextra -pedantic
    export CILKFLAG   ?= -fcilkplus -lcilkrts
    export OMPFLAG    ?= -fopenmp

    ifeq ($(TARGETARCH),POWER)
      export HTMFLAG  ?= -mhtm
      export CPUARCH ?= -mcpu=native
    else
      export HTMFLAG  ?= -mrtm
      export CPUARCH ?= -march=native
    endif
  endif
endif

ifeq ($(TOOLSET),)
  ifneq (,$(findstring clang,$(COMPILERVERSION)))
    export TOOLSET    := CLANG
    export OPTFLAG    ?= -O2
    export DBGOPTFLAG ?= -O0
    export THREADFLAG ?= -pthread
    export WARNFLAG   ?= -Wall -Wextra -pedantic
    export OMPFLAG    ?= -fopenmp

    ifeq ($(TARGETARCH),POWER)
      export HTMFLAG  ?= -mhtm
      export CPUARCH ?= -mcpu=native
    else
      export HTMFLAG  ?= -mrtm
      export CPUARCH ?= -march=native
    endif
  endif
endif

ifeq ($(TOOLSET),)
  ifneq (,$(findstring Intel,$(COMPILERVERSION)))
    export TOOLSET    := ICC
    export OPTFLAG    ?= -O2
    export DBGOPTFLAG ?= -O0
    export CPUARCH    ?= -march=native
    export THREADFLAG ?= -pthread
    export WARNFLAG   ?= -Wall -Wextra -pedantic
    export HTMFLAG    ?= -mrtm
    export CILKFLAG   ?= -lcilkrts
    export OMPFLAG    ?= -fopenmp
  endif
endif

ifeq ($(TOOLSET),)
  ifneq (,$(findstring IBM,$(COMPILERVERSION)))
    export TOOLSET    := XLC
    export OPTFLAG    ?= -O3
    export DBGOPTFLAG ?=
    export CPUARCH    ?= -qarch=auto
    export THREADFLAG ?= -pthread
    export OMPFLAG    ?= -fopenmp
  endif
endif

ifeq ($(TOOLSET),)
  ifneq (,$(findstring Sun,$(COMPILERVERSION)))
    export TOOLSET    := SUN
    export OPTFLAG    ?= -fast
    export DBGOPTFLAG ?=
    export CPUARCH    ?= -native
    export THREADFLAG ?= -mt
    export WARNFLAG   ?= -Wall -Wextra -pedantic
    export OMPFLAG    ?= -xopenmp=parallel
  endif
endif

ifeq ($(TOOLSET),)
  ifneq (,$(findstring PGI,$(COMPILERVERSION)))
    export TOOLSET    := PGI
    export OPTFLAG    ?= -O2
    export CPUARCH    ?= -ta=multicore
    export OMPFLAG    ?= -mp
    $(error PGI is currently not supported)
  endif
endif


$(info *** UCL_HOME set to default: $(UCL_HOME))
$(info *** toolset: $(TOOLSET))
$(info *** target system: $(TARGETARCH))

ifeq ($(DBGFLAG),)
  DBGFLAG:=-DNDEBUG=1
else
  OPTFLAG:=$(DBGOPTFLAG)
  $(info *** debug-mode: changed OPTFLAG to $(DBGOPTFLAG))
endif
