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

set -e -x
CMESS=$(pwd)

yum clean all
yum makecache fast
yum install -y cmake3 bzip2-devel zlib-devel xz-devel arpack-devel suitesparse-devel lapack-devel blas-devel matio-devel hdf5-devel SuperLU-devel
yum install -y scipy numpy
cd ..
mkdir -p build
cd build
cmake3 ${CMESS} -DDEBUG=OFF -DARPACK=ON -DMATIO=ON -DPYTHON=ON -DOPENMP=OFF
make
python  python/setup.2.6.py build && python  python/setup.2.6.py install
cd tests/matrix
make test

