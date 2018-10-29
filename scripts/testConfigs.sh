#!/bin/bash

test_make()
{
  selector=$1
  container=$2
  allocator=$3
  expected=$4

  outdir=${COMPDIR:-"./tmp"}

  echo "make TEST_CONTAINER=$container TEST_ALLOC=$allocator -f Makefile.$selector"
  make TEST_CONTAINER="$container" TEST_ALLOC="$allocator" -f Makefile."$selector"
  res="$?"

  if [[ "$expected" -ne "$res" ]]; then
    exit $res
  fi

  if [[ "$expected" -eq 0 ]]; then
    echo run
    mv "$outdir"/"$selector".bin "$outdir"/"$selector""-$2-$3-$CXX.bin"
  fi
}

test_simple_make()
{
  selector=$1
  testname=$2
  expected=$3

  outdir=${OUTPUTDIR:-"./tmp"}

  echo "make TEST_NAME=$testname -f Makefile.$1"
  make "TEST_NAME=$testname" -f Makefile.$1
  res="$?"

  if [[ "$expected" -ne "$res" ]]; then
    exit $res
  fi

  if [[ "$expected" -eq 0 ]]; then
    echo run
    mv "$outdir"/"$selector""Tester.bin" "$outdir"/"$selector""Tester-""$testname""-$CXX.bin"
  fi
}


test_htm_make()
{
  if [[ $host == "iprogress-phi" ]] && [[ $CXX == "icpc" ]]; then
    return
  fi

  if [[ $host == "blueblaze" ]] && [[ $CXX == "xlc++" ]]; then
    return
  fi

  selector=$1
  allocator=$2
  expected=$3

  outdir=${COMPDIR:-"./tmp"}

  echo "make TEST_ALLOC=$allocator TEST_HTM=1 -f Makefile.$selector"
  make TEST_ALLOC="$allocator" TEST_HTM=1 -f Makefile."$selector"
  res="$?"

  if [[ "$expected" -ne "$res" ]]; then
    exit $res
  fi

  if [[ "$expected" -eq 0 ]]; then
    echo run
    mv "$outdir"/"$selector".bin "$outdir"/"$selector""-HTM-$2-""$CXX"".bin"
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
  test_simple_make stack TEST_NO_MANAGER 0
  test_simple_make stack TEST_EPOCH_MANAGER 0
  test_simple_make stack TEST_PUB_SCAN_MANAGER 0
  test_simple_make stack TEST_GC_MANAGER 0
  test_simple_make stack TEST_STD_LOCKGUARD 0
  test_simple_make stack TEST_UAB_LOCKGUARD 0
  #~ test_stack_make TEST_UAB_ELIDEGUARD 0
}

test_queues()
{
  test_simple_make queue TEST_NO_MANAGER 0
  test_simple_make queue TEST_EPOCH_MANAGER 0
  test_simple_make queue TEST_PUB_SCAN_MANAGER 0
  test_simple_make queue TEST_GC_MANAGER 0
  test_simple_make queue TEST_STD_LOCKGUARD 0
  test_simple_make queue TEST_UAB_LOCKGUARD 0
  #~ test_stack_make TEST_UAB_ELIDEGUARD 0
}



###
# COMPILERS
COMPILERS="$COMPILERS g++ g++-5 g++-6 g++-7 g++-8"
COMPILERS="$COMPILERS clang++ clang++-3.4 clang++-3.5 clang++-3.6 clang++-3.7 clang++-4.0 clang++-5.0 clang++-6.0"
COMPILERS="$COMPILERS icpc xlc++ sunCC"

###
# DATA STRUCTURE TESTS

# TESTS="test_queues"

TESTS="test_skiplists test_stacks test_queues"

for arg in $COMPILERS
do
  if hash $arg 2>/dev/null; then
    export CXX=$arg
    echo "* testing $CXX"
    for test in $TESTS
    do
      $test
    done
  fi
done
