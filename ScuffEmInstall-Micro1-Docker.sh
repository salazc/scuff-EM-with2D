#!/bin/bash


# Cleaning Repo
git clean -dfx

apt-get install gfortran build-essential g++ gcc
apt-get install libblas-dev liblapack-dev libgsl-dev

apt-get install libopenblas-dev libhdf5-openmpi-dev python3-scipy gmsh
apt-get install python-dev-is-python3

apt-get install autoconf
apt-get install swig
apt-get install flex bison cmake make

apt-get install libtool autoconf automake pkg-config libtinyxml2-dev

apt-get install python3-pip python3-matplotlib python3-numpy python3-opencv

apt-get update

# PyTorch (CPU) is not in Ubuntu 24.04 repos; install via pip
pip3 install --no-cache-dir --break-system-packages \
    torch --index-url https://download.pytorch.org/whl/cpu

export CC=mpicc
export CXX=mpic++

# This may not be needed but it generates the parser.c file
cd libs/libMatProp/cmatheval
./RunLexAndYacc.sh
cd ../../../



sh ./autogen.sh  --prefix=${HOME}/scuff-em-installation  CPPFLAGS="-I/usr/include/hdf5/openmpi" LDFLAGS="-L/usr/lib/x86_64-linux-gnu/hdf5/openmpi -L/usr/lib/x86_64-linux-gnu" LIBS="-ltinyxml2"

make -j 8 install

