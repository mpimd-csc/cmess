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
apt-get install --yes cmake3 libopenblas-dev libsuitesparse-dev zlib1g-dev libbz2-dev liblzma-dev
apt-get install --yes libarpack2-dev libmatio-dev libsuperlu3-dev liblapack-dev
apt-get install --yes python-dev python3-dev python-pip python3-pip python-setuptools python3-setuptools
pip2 install --upgrade pip numpy scipy
pip3 install --upgrade pip numpy scipy coverage

cd ..
mkdir -p build
cd build
cmake ${CMESS} -DDEBUG=OFF -DARPACK=ON -DMATIO=ON -DPYTHON=ON -DOPENMP=OFF
make

# build, install and test python
PYTHON_VERSION=$(python -c "print '%d.%d'%( __import__('sys').version_info[:2])")
PYTHON3_VERSION=$(python3 -c "print('%d.%d'%( __import__('sys').version_info[:2]))")

make pymess-${PYTHON_VERSION}-build
make pymess-${PYTHON_VERSION}-install
make pymess-${PYTHON_VERSION}-run-examples
make pymess-${PYTHON_VERSION}-run-tests

make pymess-${PYTHON3_VERSION}-build
make pymess-${PYTHON3_VERSION}-install
make pymess-${PYTHON3_VERSION}-run-examples
make pymess-${PYTHON3_VERSION}-run-tests

# run python 3 with coverage
make pymess-${PYTHON3_VERSION}-coverage


