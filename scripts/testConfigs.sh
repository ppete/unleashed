#!/bin/bash

test_make()
{
  selector=$1
  container=$2
  allocator=$3
  expected=$4

  make TEST_CONTAINER="$container" TEST_ALLOC="$allocator" -f Makefile."$selector"
  res="$?"

  if [[ "$expected" -ne "$res" ]]; then
    exit $res
  fi

  if [[ "$expected" -eq 0 ]]; then
    echo run
    mv ./tmp/"$selector".bin ./tmp/"$selector""-$2-$3-$CXX.bin"
  fi
}

test_htm_make()
{
  selector=$1
  allocator=$2
  expected=$3

  make TEST_ALLOC="$allocator" TEST_HTM=1 -f Makefile."$selector"
  res="$?"

  if [[ "$expected" -ne "$res" ]]; then
    exit $res
  fi

  if [[ "$expected" -eq 0 ]]; then
    echo run
    mv ./tmp/"$selector".bin ./tmp/"$selector""-HTM-$2-""$CXX"".bin"
  fi
}


test_skiplists()
{
  test_make skiplist TEST_LOCKFREE_SKIPLIST TEST_NO_MANAGER 0
  test_make skiplist TEST_LOCKFREE_SKIPLIST TEST_EPOCH_MANAGER 0
  test_make skiplist TEST_LOCKFREE_SKIPLIST TEST_PUB_SCAN_MANAGER 0
  test_make skiplist TEST_LOCKFREE_SKIPLIST TEST_GC_MANAGER 0

  test_make skiplist TEST_LOCKING_SKIPLIST TEST_NO_MANAGER 0
  test_make skiplist TEST_LOCKING_SKIPLIST TEST_EPOCH_MANAGER 0
  test_make skiplist TEST_LOCKING_SKIPLIST TEST_PUB_SCAN_MANAGER 2
  test_make skiplist TEST_LOCKING_SKIPLIST TEST_GC_MANAGER 0

  test_htm_make skiplist TEST_NO_MANAGER 0
  test_htm_make skiplist TEST_EPOCH_MANAGER 0
  test_htm_make skiplist TEST_REFCOUNT_MANAGER 0
  test_htm_make skiplist TEST_PUB_SCAN_MANAGER 0
  test_htm_make skiplist TEST_STACKTRACK_MANAGER 0
}

test_stacks()
{
  test_make stack TEST_LOCKFREE_STACK TEST_NO_MANAGER 0
  test_make stack TEST_LOCKFREE_STACK TEST_EPOCH_MANAGER 0
  test_make stack TEST_LOCKFREE_STACK TEST_PUB_SCAN_MANAGER 0
  test_make stack TEST_LOCKFREE_STACK TEST_GC_MANAGER 0

  #~ test_make stack TEST_LOCKING_STACK TEST_DEFAULT_ALLOC 0
  #~ test_htm_make stack TEST_NO_MANAGER 0
  #~ test_htm_make stack TEST_EPOCH_MANAGER 0
  #~ test_htm_make stack TEST_REFCOUNT_MANAGER 0
  #~ test_htm_make stack TEST_PUB_SCAN_MANAGER 0
  #~ test_htm_make stack TEST_STACKTRACK_MANAGER 0
}


###
# COMPILERS
COMPILERS=

# GNU compilers
# COMPILERS="g++-4.8 g++-4.9 g++-5.0 powerpc-linux-gnu-g++-4.9 arm-linux-gnueabihf-g++-4.9 x86_64-w64-mingw32-g++"
COMPILERS="$COMPILERS g++"

# Clang family
# COMPILERS="$COMPILERS clang++-3.4 clang++-3.5 clang++-3.6 clang++-3.7"

# Intel compiler
# COMPILERS="$COMPILERS icpc"

# IBM compilers
# 
# COMPILERS="$COMPILERS xlc++"


###
# DATA STRUCTURE TESTS

#~ TESTS="test_skiplists test_stacks"

TESTS="test_skiplists"

for arg in $COMPILERS
do
  if hash $arg 2>/dev/null; then
    export CXX=$arg
    for test in $TESTS
    do
      $test
    done
  fi
done
