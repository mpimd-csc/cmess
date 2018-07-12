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
apt-get install --yes libarpack2-dev libmatio-dev libsuperlu3-dev
apt-get install --yes python-dev python-scipy python-numpy
apt-get install --yes python3-dev python3-scipy python3-numpy
cd ..
mkdir -p build
cd build
cmake ${CMESS} -DDEBUG=OFF -DARPACK=ON -DMATIO=ON -DPYTHON=ON -DOPENMP=OFF
make
export LD_LIBRARY_PATH=$(pwd)/lib/:$LD_LIBRARY_PATH
python  python/setup.2.7.py build && python  python/setup.2.7.py install
python3 python/setup.3.2.py build && python3 python/setup.3.2.py install

#python tests exclude CALLBACK because scipy is too old from apt
python  python/unittests/run.py EASY
python  python/unittests/run.py CONVERSION
python  python/unittests/run.py LRADI
python  python/unittests/run.py LRNM


#python3 tests exclude CALLBACK because scipy is too old from apt
python3  python/unittests/run.py EASY
python3  python/unittests/run.py CONVERSION
python3  python/unittests/run.py LRADI
python3  python/unittests/run.py LRNM


# run examples
#PYTHON_VERSION=$(python -c "print '%d.%d.%d'%( __import__('sys').version_info[:3])")
#cd python/example
#./run_${PYTHON_VERSION}.sh
