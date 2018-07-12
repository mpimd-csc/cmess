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
#apt-get install --yes python-pip
#pip install --upgrade pip
easy_install pip
pip install pylint

# generate build
cd ..
mkdir -p build
cd build
cmake ${CMESS} -DDEBUG=OFF -DARPACK=ON -DMATIO=ON -DPYTHON=ON -DOPENMP=OFF

# compile cmess and pymess
make

# install pymess
make pymess-2.7-build
make pymess-2.7-install

# run pylint on pymess
make pymess-lint
