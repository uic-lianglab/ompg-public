#!/bin/bash

###################################################
# Script paths
###################################################

# Path to script directory
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"

###################################################
# Build NAMD
###################################################

# Store current directory
pushd .

# Switch to script directory
cd $SCRIPT_DIR

# Unpack matching Charm++ source code:
tar xf charm-6.7.0.tar

# Build and test the Charm++/Converse library (single-node multicore version):
cd charm-6.7.0
./build charm++ multicore-linux64 --with-production
cd multicore-linux64/tests/charm++/megatest
make pgm
./pgm +p4

# Install TCL and FFTW libraries:
cd $SCRIPT_DIR
tar xzf fftw-linux-x86_64.tar.gz
mv linux-x86_64 fftw
tar xzf tcl8.5.9-linux-x86_64.tar.gz
tar xzf tcl8.5.9-linux-x86_64-threaded.tar.gz
mv tcl8.5.9-linux-x86_64 tcl
mv tcl8.5.9-linux-x86_64-threaded tcl-threaded

# Set up build directory and compile:
./config Linux-x86_64-g++ --charm-arch multicore-linux64
cd Linux-x86_64-g++
make

# Restore directory
popd

echo Finished

