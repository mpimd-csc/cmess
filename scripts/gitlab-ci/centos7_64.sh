#!/bin/bash

#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>.
#
# Copyright (C) Peter Benner, Martin Koehler, Jens Saak and others
#               2009-2018
#

set -e
CMESS=$(pwd)

yum makecache fast
#install 64 bit blas lapack
yum install -y lapack64-devel.x86_64 blas64-devel.x86_64

#install other packages
yum install -y cmake3 bzip2-devel zlib-devel xz-devel  hdf5-devel wget

#download, compile, link against 64 library, install suitesparse and delete sources
CMESS=$(pwd)
cd ..
wget http://faculty.cse.tamu.edu/davis/SuiteSparse/SuiteSparse-4.5.3.tar.gz
tar xfv SuiteSparse-4.5.3.tar.gz
mv SuiteSparse SuiteSparseSRC
cd SuiteSparseSRC
make BLAS=-lblas64 LAPACK=-llapack64 OPTIMIZATION=-g  CF="-g -fexceptions -fPIC -fopenmp -DLONGBLAS=long" F77FLAGS+="-fdefault-integer-8" MKLROOT= INSTALL=/SuiteSparse/ CHOLMOD_CONFIG=-DNPARTITION install
cd ..
rm -rf SuiteSparseSRC SuiteSparse-4.5.3.tar.gz

#perform cmess 64 bit build
cd ${CMESS}
cd ..
mkdir -p build
cd build
cmake3 ${CMESS} -DDEBUG=OFF -DMESS64=ON -DBLAS=/usr/lib64/libblas64.so -DLAPACK=/usr/lib64/liblapack64.so -DSUITESPARSE=/SuiteSparse/
make
make test
