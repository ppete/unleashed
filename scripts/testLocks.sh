#!/usr/bin/env bash

test_make()
{
  lock=$1
  expected=$2

  echo "make TARGET=$COMPDIR/test-$lock-$CXX-$CXXVERSION.bin TEST_LOCK=$lock -f Makefile.spinlock -C ../examples/locks"
  make TARGET="$COMPDIR"/test-"$lock"-"$CXX"-"$CXXVERSION".bin TEST_LOCK="$lock" -f Makefile.spinlock -C ../examples/locks
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

# GNU compilers
# COMPILERS="g++-4.8 g++-4.9 g++-5.0 powerpc-linux-gnu-g++-4.9 arm-linux-gnueabihf-g++-4.9 x86_64-w64-mingw32-g++"

# Clang family
# COMPILERS="$COMPILERS clang++-3.4 clang++-3.5 clang++-3.6 clang++-3.7"

# Intel compiler
# COMPILERS="$COMPILERS icpc"


###
# COMPILERS
#~ COMPILERS="$COMPILERS g++-5 g++-6 g++-7 g++-8"
#~ COMPILERS="$COMPILERS clang++-4.0 clang++-6.0 clang++-7"
#~ COMPILERS="$COMPILERS icpc xlc++ sunCC"
if [ -z ${COMPILERS+x} ]; 
then
  COMPILERS="g++ clang++ icpx"
fi

STANDARDS="-std=c++11 -std=c++20 -std=c++14 -std=c++17"

###
# DATA STRUCTURE TESTS

#~ TESTS="test_queues"

#~ TESTS="test_skiplists test_stacks test_queues"
TESTS="test_locks"

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
