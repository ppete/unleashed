#!/bin/bash

test_make()
{
  selector=$1
  container=$2
  allocator=$3
  expected=$4

  outdir=$COMPDIR

  echo "make TEST_CONTAINER=$container TEST_ALLOC=$allocator -f Makefile.$selector"
  make TEST_CONTAINER="$container" TEST_ALLOC="$allocator" -f Makefile."$selector" -C ../examples/containers
  res="$?"

  if [[ "$expected" -ne "$res" ]]; then
    exit $res
  fi

  if [[ "$expected" -eq 0 ]]; then
    mv "$outdir"/"$selector".bin "$outdir"/"$selector""-$2-$3-$CXX.bin"
  fi
}

test_simple_make()
{
  selector=$1
  testname=$2
  expected=$3

  outdir="$COMPDIR"

  echo "make TEST_NAME=$testname TARGET=$outdir/$selector-$2-$3-$comp-$CXXVERSION.bin -f Makefile.$selector -C ../examples/containers"
  make TEST_NAME="$testname" TARGET="$outdir/$selector-$2-$3-$comp-$CXXVERSION.bin" -f Makefile."$selector" -C ../examples/containers
  res="$?"

  if [[ "$expected" -ne "$res" ]]; then
    exit $res
  fi
  
#  if [[ "$expected" -eq 0 ]]; then
#    mv "$outdir"/"$selector".bin "$outdir"/"$selector""-$2-$3-$CXX-$CXXVERSION.bin"
#  fi
}


test_htm_make()
{
  selector=$1
  allocator=$2
  expected=$3

  outdir="$COMPDIR"

  echo "make TEST_ALLOC=$allocator TEST_HTM=1 -f Makefile.$selector"
  make TEST_ALLOC="$allocator" TEST_HTM=1 -f Makefile."$selector"
  res="$?"

  if [[ "$expected" -ne "$res" ]]; then
    exit $res
  fi

  if [[ "$expected" -eq 0 ]]; then
    mv "$outdir"/"$selector".bin "$outdir"/"$selector""-HTM-$2-""$CXX"".bin"
  fi
}


test_skiplists()
{
  test_make skiplist TEST_LF_CONTAINER TEST_NO_MANAGER 0
  test_make skiplist TEST_LF_CONTAINER TEST_EPOCH_MANAGER 0
  test_make skiplist TEST_LF_CONTAINER TEST_PUB_SCAN_MANAGER 0
  test_make skiplist TEST_LF_CONTAINER TEST_GC_MANAGER 0

  test_make skiplist TEST_LOCK_CONTAINER TEST_NO_MANAGER 0
  test_make skiplist TEST_LOCK_CONTAINER TEST_EPOCH_MANAGER 0
  test_make skiplist TEST_LOCK_CONTAINER TEST_PUB_SCAN_MANAGER 2
  test_make skiplist TEST_LOCK_CONTAINER TEST_GC_MANAGER 0

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
  #~ test_simple_make stack TEST_GC_MANAGER 0
  test_simple_make stack TEST_STD_LOCKGUARD 0
  #~ test_simple_make stack TEST_UNLEASHED_LOCKGUARD 0
  #~ test_stack_make TEST_UNLEASHED_ELIDEGUARD 0
}

test_queues()
{
  test_simple_make queue TEST_NO_MANAGER 0
  test_simple_make queue TEST_EPOCH_MANAGER 0
  test_simple_make queue TEST_PUB_SCAN_MANAGER 0
  #~ test_simple_make queue TEST_GC_MANAGER 0
  test_simple_make queue TEST_STD_LOCKGUARD 0
  #~ test_simple_make queue TEST_UNLEASHED_LOCKGUARD 0
  #~ test_stack_make TEST_UNLEASHED_ELIDEGUARD 0
}



###
# COMPILERS
#~ COMPILERS="$COMPILERS g++-5 g++-6 g++-7 g++-8"
#~ COMPILERS="$COMPILERS clang++-4.0 clang++-6.0 clang++-7"
#~ COMPILERS="$COMPILERS icpc xlc++ sunCC"
if [ -z ${COMPILERS+x} ]; 
then
  COMPILERS="g++-12 clang++-15 icpx"
fi

STANDARDS="-std=c++11 -std=c++14 -std=c++17 -std=c++20"

###
# DATA STRUCTURE TESTS

#~ TESTS="test_queues"

#~ TESTS="test_skiplists test_stacks test_queues"
TESTS="test_stacks test_queues"

for arg in $COMPILERS
do
  if hash $arg 2>/dev/null; then
    for std in $STANDARDS
    do
      export CXX="$arg"
      export CXXVERSION="$std"
      echo "* testing $CXX"
      for test in $TESTS
      do
        $test
      done
    done
  fi
done
