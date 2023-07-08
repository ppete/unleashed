
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

  ifneq (,$(findstring mips64,$(THISSYSTEM)))
    TARGETARCH:=MIPS64
  endif

  ifneq (,$(findstring riscv64,$(THISSYSTEM)))
    TARGETARCH:=RISCV64
  endif

  ifneq (,$(findstring arm64,$(THISSYSTEM)))
    TARGETARCH:=ARM64
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

    TESTNATIVE := 1

    ifeq ($(TARGETARCH),POWER)
      ARCHFLAG := -mcpu
    else ifeq ($(TARGETARCH),SPARC)
      ARCHFLAG := -mcpu
    else
      # x86, arm
      ARCHFLAG := -march

      export HTMFLAG  ?= -mrtm
    endif
  endif
endif

ifeq ($(TOOLSET),)
  ifneq (,$(findstring Intel,$(COMPILERVERSION)))
    export TOOLSET    := ICC
    export OPTFLAG    ?= -O2
    export DBGOPTFLAG ?= -O0
    export THREADFLAG ?= -pthread
    export WARNFLAG   ?= -Wall -Wextra -pedantic
    export HTMFLAG    ?= -mrtm
    export CILKFLAG   ?= -lcilkrts
    export OMPFLAG    ?= -fopenmp

    TESTNATIVE := 1
    ARCHFLAG := -march
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
    export CILKFLAG   ?= -fopencilk

    TESTNATIVE := 1

    ifeq ($(TARGETARCH),POWER)
      export HTMFLAG  ?= -mhtm

      ARCHFLAG := -mcpu
    else
      # Clang10 does not accept -mcpu=native on Sparc and MIPS
      ifneq ($(TARGETARCH),SPARC)
      ifneq ($(TARGETARCH),MIPS64)
        ARCHFLAG := -march
        export HTMFLAG  ?= -mrtm
      endif
      endif
    endif
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
    export OPTFLAG    ?= -O3
    export DBGOTPFLAG ?= -O0
    export CPUARCH    ?= -ta=multicore
    export OMPFLAG    ?= -mp
    export WARNFLAG   ?= -Wall -pedantic
    export THREADFLAG ?=
    $(info PGI may not be fully supported)
  endif
endif

##
## test for native flag support

ifeq ($(TESTNATIVE),1)
ifeq ($(CPUARCH),)
  ARCHFLAG:=$(ARCHFLAG)=native

  SUCCESS:=$(shell $(CXX) $(ARCHFLAG) -c $(UCL_HOME)/examples/make-util/test-hello.cc -o /dev/null >/dev/null 2>&1; echo $$?)

  ifeq ($(SUCCESS),0)
    export CPUARCH:=$(ARCHFLAG)
  endif
endif
endif

#
# test for C++-17 support
ifeq ($(CXX17_ENABLED),)
  ARCHFLAG:=$(ARCHFLAG)=native

  SUCCESS:=$(shell $(CXX) -std=c++17 -c $(UCL_HOME)/examples/make-util/test-hello.cc -o /dev/null >/dev/null 2>&1; echo $$?)

  ifeq ($(SUCCESS),0)
    export CXX17_ENABLED:=1
  endif
endif


##
## test for libatomic

SUCCESS:=$(shell $(CXX) -latomic $(UCL_HOME)/examples/make-util/test-hello.cc -o /dev/null >/dev/null 2>&1; echo $$?)

ifeq ($(SUCCESS),0)
  export LIBATOMIC:=-latomic
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
