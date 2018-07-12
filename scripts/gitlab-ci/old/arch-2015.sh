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

pacman -Sy
pacman -S --noconfirm suitesparse blas lapack libmatio arpack zlib bzip2
pacman -S --noconfirm python2 python2-numpy python2-scipy
pacman -S --noconfirm python python-numpy python-scipy
cd ..
mkdir -p build
cd build
cmake ../cmess -DDEBUG=OFF -DARPACK=ON -DMATIO=ON -DPYTHON=ON -DOPENMP=OFF
make
export LD_LIBRARY_PATH=$(pwd)/lib/:$LD_LIBRARY_PATH
python2 python/setup.2.7.py build && python2 python/setup.2.7.py install
python  python/setup.3.5.py build && python  python/setup.3.5.py install
cd tests/matrix
make test

