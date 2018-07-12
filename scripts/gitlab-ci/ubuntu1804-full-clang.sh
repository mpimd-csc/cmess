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

apt-get update --yes
apt-get install --yes cmake libopenblas-dev libsuitesparse-dev zlib1g-dev libbz2-dev liblzma-dev
apt-get install --yes libarpack2-dev libmatio-dev libsuperlu-dev


# install different clang versions
apt-get install --yes   wget software-properties-common
wget -O - https://apt.llvm.org/llvm-snapshot.gpg.key| apt-key add -
apt-add-repository "deb http://apt.llvm.org/bionic/ llvm-toolchain-bionic main"
apt-add-repository "deb http://apt.llvm.org/bionic/ llvm-toolchain-bionic-5.0 main"
apt-add-repository "deb http://apt.llvm.org/bionic/ llvm-toolchain-bionic-6.0 main"
apt-get update
apt-get install --yes clang-3.9 clang-4.0 clang-5.0 clang-6.0



# compile for each version
for VERSION in  3.9 4.0 5.0 6.0
do
    echo "COMPILE WITH CLANG VERSION: ${VERSION}"                                   &&  \
    cd ..                                                                           &&  \
    mkdir -p build${VERSION}                                                        &&  \
    cd build${VERSION}                                                              &&  \
    export CC=clang-${VERSION} CXX=clang++-${VERSION}                               &&  \
    cmake ${CMESS} -DDEBUG=OFF -DARPACK=ON -DMATIO=ON -DPYTHON=OFF -DOPENMP=OFF     &&  \
    make                                                                            &&  \
    cd ..
done



