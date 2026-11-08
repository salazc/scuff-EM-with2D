
# Cleaning Repo
git clean -dfx

sudo apt install gfortran build-essential g++ gcc
sudo apt install libblas-dev liblapack-dev libgsl-dev

sudo apt-get install libopenblas-dev libhdf5-openmpi-dev python3-scipy gmsh
sudo apt-get install python-dev-is-python3

sudo apt-get install autoconf
sudo apt-get install swig
sudo apt-get install flex bison cmake make

sudo apt install libtool autoconf automake pkg-config libtinyxml2-dev

apt-get update

export CC=mpicc
export CXX=mpic++

# This may not be needed but it generates the parser.c file
cd libs/libMatProp/cmatheval
./RunLexAndYacc.sh
cd ../../../



sh ./autogen.sh  --prefix=/scuff-em-installation  CPPFLAGS="-I/usr/include/hdf5/openmpi" LDFLAGS="-L/usr/lib/x86_64-linux-gnu/hdf5/openmpi -L/usr/lib/x86_64-linux-gnu" LIBS="-ltinyxml2"

make -j 8 install

