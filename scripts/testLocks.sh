#!/usr/bin/env bash

test_make()
{
  lock=$1
  expected=$2

  echo "$MAKE TARGET=$COMPDIR/test-$lock-$CXX-$CXXVERSION.bin TEST_LOCK=$lock -f Makefile.spinlock -C ../examples/locks"
  $MAKE TARGET="$COMPDIR"/test-"$lock"-"$CXX"-"$CXXVERSION".bin TEST_LOCK="$lock" -f Makefile.spinlock -C ../examples/locks
  res="$?"

  if [[ "$expected" -ne "$res" ]]; then
    exit $res
  fi
}

test_locks()
{
  test_make TEST_ANDERSON 0
  test_make TEST_COUNTING 0
  test_make TEST_MCS 0
  test_make TEST_CLH 0
  test_make TEST_TTAS 0
  test_make TEST_TTAS_BO 0
  test_make TEST_MUTEX 0
}


###
# COMPILERS

TESTS="test_locks"

if [ -z ${MAKE+x} ]; 
then
  MAKE="make"
fi


if [ -z ${COMPILERS+x} ]; 
then
  COMPILERS="g++ clang++ icpx"
fi

if [ -z ${STANDARDS+x} ]; 
then
  STANDARDS="-std=c++11 -std=c++20 -std=c++14 -std=c++17"
fi

###
# DATA STRUCTURE TESTS

#~ TESTS="test_queues"

#~ TESTS="test_skiplists test_stacks test_queues"

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
