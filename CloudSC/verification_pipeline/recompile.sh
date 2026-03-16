#!/bin/bash
set -ex
c++ -O3 -g -std=c++20 -DNDEBUG  -Wall -Wextra -Wno-parentheses-equality -Wno-unused-variable -Wno-unused-label -Icodegen -Iinclude -I/home/primrose/Work/dace/dace/runtime/include -I/usr/include/hdf5/serial cloudsc_main.cpp codegen/*.cpp -o cloudsc_cpu_bin -lpthread -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5 || exit 1
