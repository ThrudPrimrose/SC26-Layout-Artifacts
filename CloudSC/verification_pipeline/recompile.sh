#!/bin/bash
set -e
c++ -O3 -std=c++20 -DNDEBUG -Wno-parentheses-equality \
  -Icodegen -Iinclude \
  -I$HOME/Work/dace/dace/runtime/include \
  -I/usr/include/hdf5/serial \
  cloudsc_main.cpp codegen/*.cpp -o cloudsc_cpu_bin \
  -lpthread \
  -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5