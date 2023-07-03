# Tests and Benchmarks for UCL Tasks

The tests are available for a variety of task systems (i.e., UCL, OpenMP, TBB, OpenCilk). 

Some of the tests were adapted from the Barcelona OpenMP Task Test Suite (BOTS). The BOTS tests are distributed under a *GNU GPL license* - see readme files in the corresponding test directories.

Compile the tests with make (or gmake on Unix/BSD platforms).

```bash
make
```

The make command supports the following flags:
- CXX: set the compiler (e.g., g++-12, clang++, icpx, ..). The Makefile will run some tests to figure out which programming models are supported and what default flags to use. The environment variables (TBB_HOME and QTHREADS_HOME) need to point to the installation directory for Threading Building Blocks and QThreads respectively.
- OPTFLAG: customize optimization flag
- ARCHFLAG: customize architecture flag
- WARNFLAG: customize warning level

```bash
make CXX=clang++-15 OPTFLAG=-O3 
```


