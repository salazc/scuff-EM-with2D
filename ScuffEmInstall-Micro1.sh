

apt-get install libopenblas-dev libhdf5-openmpi-dev python-dev python3-scipy gmsh
apt-get install python-dev-is-python3

apt-get install autoconf
apt-get install swig
apt-get install flex bison

apt-get update

export CC=mpicc
export CXX=mpic++

# This may not be needed but it generates the parser.c file
cd libs/libMatProp/cmatheval
./RunLexAndYacc.sh
cd ../../../



sh ./autogen.sh  --prefix=/scuff-em-installation  CPPFLAGS=-I/usr/include/hdf5/openmpi LDFLAGS=-L/usr/lib/x86_64-linux-gnu/hdf5/openmpi LIBS=-ltinyxml2


make -j 8 install

