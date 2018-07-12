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

# stop on error
set -e -x
CMESS=$(pwd)

#some links to packages
NUMPY_PY27=http://pc1184.mpi-magdeburg.mpg.de:8080/soft/numpy-1.13.3-cp27-cp27mu-manylinux1_x86_64.whl
NUMPY_PY34=http://pc1184.mpi-magdeburg.mpg.de:8080/soft/numpy-1.13.3-cp34-cp34m-manylinux1_x86_64.whl
NUMPY_PY35=http://pc1184.mpi-magdeburg.mpg.de:8080/soft/numpy-1.13.3-cp35-cp35m-manylinux1_x86_64.whl

SCIPY_PY27=http://pc1184.mpi-magdeburg.mpg.de:8080/soft/scipy-0.16.0-cp27-cp27mu-manylinux1_x86_64.whl
SCIPY_PY34=http://pc1184.mpi-magdeburg.mpg.de:8080/soft/scipy-0.16.0-cp34-cp34m-manylinux1_x86_64.whl
SCIPY_PY35=http://pc1184.mpi-magdeburg.mpg.de:8080/soft/scipy-0.16.0-cp35-cp35m-manylinux1_x86_64.whl


#check if MESS_ARTIFACTS variable is defined from .gitlab-ci.yml
if [ -z "${MESS_ARTIFACTS}" ]
    then
    echo  "MESS_ARTIFACTS variable is not defined"
    DEPLOY=${CMESS}/mess_artifacts
    else
    echo  "MESS_ARTIFACTS variable is defined"
    DEPLOY=${CMESS}/${MESS_ARTIFACTS}
fi


apt-get update --yes
apt-get install --yes cmake libopenblas-dev zlib1g-dev libbz2-dev liblzma-dev libsuitesparse-dev
apt-get install --yes python-dev python3-dev python-setuptools python3-setuptools
easy_install pip
easy_install3 pip


# update pip, install numpy scipy
pip2 install --upgrade pip
pip2 install ${NUMPY_PY27} ${SCIPY_PY27} wheel
pip3 install --upgrade pip
pip3 install ${NUMPY_PY35} ${SCIPY_PY35} wheel


# install the pymess wheels
cd ${DEPLOY}
pip2 install *cp27*.whl
pip3 install *cp35*.whl


# configure C-M.E.S.S. build and NOT compile or install anything
cd ${CMESS}
cd ..
mkdir build
cd build
cmake ${CMESS} -DPYTHON=ON

# python examples are configured, we can test our python wheels now
make pymess-2.7-run-examples
make pymess-3.5-run-examples




