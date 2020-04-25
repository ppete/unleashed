#!/usr/bin/env bash

test_make()
{
  lock=$1
  expected=$2

  make TEST_LOCK="$lock" -f Makefile.spinlock
  res="$?"

  if [[ "$expected" -ne "$res" ]]; then
    exit $res
  fi
}

test_locks()
{
  cd ../examples/locks

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
# LOCK TESTS
test_locks

