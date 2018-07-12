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
yum install -y cmake3 bzip2-devel zlib-devel xz-devel arpack-devel suitesparse-devel lapack-devel blas-devel hdf5-devel SuperLU-devel
yum install -y matio-devel
yum install -y scipy numpy
cd ..
mkdir -p build
cd build
cmake ${CMESS} -DDEBUG=OFF -DARPACK=ON -DMATIO=ON -DPYTHON=ON -DOPENMP=OFF
make


PYTHON_VERSION=$(python -c "print '%d.%d'%( __import__('sys').version_info[:2])")
make pymess-${PYTHON_VERSION}-build
make pymess-${PYTHON_VERSION}-install


cd tests/matrix
make test

