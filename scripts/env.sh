#!/usr/bin/env bash

# common environment variables for tests

# MAKE points to gnu make
export MAKE=`which make`

export COMPDIR=/PATH/TO/COMPILER/OUTPUT

# popular compilers
export COMPILERS="g++ clang++"

# C++ versions
export STANDARDS="-std=c++11 -std=c++20"
