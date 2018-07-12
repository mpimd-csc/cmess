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

#check if MESS_ARTIFACTS variable is defined from .gitlab-ci.yml
if [ -z "${MESS_ARTIFACTS}" ]
    then
    echo  "MESS_ARTIFACTS variable is not defined"
    DEPLOY=${CMESS}/mess_artifacts
    else
    echo  "MESS_ARTIFACTS variable is defined"
    DEPLOY=${CMESS}/${MESS_ARTIFACTS}
fi

zypper -n install cmake libbz2-devel zlib-devel xz-devel  suitesparse-devel lapack-devel blas-devel
zypper -n install python-devel python3-devel
zypper -n install python-setuptools python3-setuptools python-pip python3-pip

pip2 install --upgrade pip
pip3 install --upgrade pip

# install the pymess wheels
cd ${DEPLOY}
pip2 install *cp27*.whl
pip3 install *cp34*.whl

cd ${CMESS}
cd ..
mkdir -p build
cd build
cmake ${CMESS} -DPYTHON=ON

PYTHON_VERSION=$(python -c "print '%d.%d'%( __import__('sys').version_info[:2])")
PYTHON3_VERSION=$(python3 -c "print('%d.%d'%( __import__('sys').version_info[:2]))")

make pymess-${PYTHON_VERSION}-run-examples
make pymess-${PYTHON3_VERSION}-run-examples


