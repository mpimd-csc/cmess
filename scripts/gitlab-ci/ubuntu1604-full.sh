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

apt-get update --yes || true
apt-get install --yes cmake libopenblas-dev libsuitesparse-dev zlib1g-dev libbz2-dev liblzma-dev
apt-get install --yes libarpack2-dev libmatio-dev libsuperlu-dev
apt-get install --yes python-dev python3-dev python-scipy python3-scipy python-setuptools python3-setuptools

# install different gcc versions
apt-get install --yes  software-properties-common  python-software-properties
add-apt-repository --yes ppa:ubuntu-toolchain-r/test
apt-get update || true
apt-get install --yes   gcc-4.9 gcc-5 gcc-6
apt-get install --yes   g++-4.9 g++-5 g++-6
apt-get install --yes   gfortran-4.9 gfortran-5 gfortran-6

# compile for each version
for VERSION in 4.9 5 6
do
    echo "COMPILE WITH GCC VERSION: ${VERSION}"                                             &&  \
    cd ..                                                                                   &&  \
    mkdir -p build${VERSION}                                                                &&  \
    cd build${VERSION}                                                                      &&  \
    export CC=gcc-${VERSION} CXX=g++-${VERSION} FC=gfortran-${VERSION}                      &&  \
                                                                                                \
    cmake ${CMESS} -DDEBUG=OFF -DARPACK=ON -DMATIO=ON -DPYTHON=ON -DOPENMP=ON               &&  \
    make                                                                                    &&  \
                                                                                                \
    PYTHON_VERSION=$(python -c "print '%d.%d'%( __import__('sys').version_info[:2])")       &&  \
    PYTHON3_VERSION=$(python3 -c "print('%d.%d'%( __import__('sys').version_info[:2]))")    &&  \
                                                                                                \
    make pymess-${PYTHON_VERSION}-build                                                     &&  \
    make pymess-${PYTHON_VERSION}-install                                                   &&  \
                                                                                                \
    make pymess-${PYTHON3_VERSION}-build                                                    &&  \
    make pymess-${PYTHON3_VERSION}-install                                                  &&  \
                                                                                                \
    cd ..
done

