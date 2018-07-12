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

SCIPY_PY27=http://pc1184.mpi-magdeburg.mpg.de:8080/soft/scipy-0.16.0-cp27-cp27mu-manylinux1_x86_64.whl
SCIPY_PY34=http://pc1184.mpi-magdeburg.mpg.de:8080/soft/scipy-0.16.0-cp34-cp34m-manylinux1_x86_64.whl


#check if MESS_ARTIFACTS variable is defined from .gitlab-ci.yml
if [ -z "${MESS_ARTIFACTS}" ]
    then
    echo  "MESS_ARTIFACTS variable is not defined"
    DEPLOY=${CMESS}/mess_artifacts
    else
    echo  "MESS_ARTIFACTS variable is defined"
    DEPLOY=${CMESS}/${MESS_ARTIFACTS}
fi


# update yum
yum clean all
yum makecache fast
yum -y update
yum install -y cmake3 bzip2-devel zlib-devel xz-devel suitesparse-devel lapack-devel blas-devel
yum install -y python2-pip python34-pip python34-devel python-devel


# install scipy numpy
pip2.7 install --upgrade pip
pip2.7 install ${NUMPY_PY27} ${SCIPY_PY27} wheel
pip3.4 install --upgrade pip
pip3.4 install ${NUMPY_PY34} ${SCIPY_PY34} wheel


# now install the deployed PYMESS
cd ${DEPLOY}
pip2.7 install *cp27*.whl
pip3.4 install *cp34*.whl


# configure C-M.E.S.S. build and NOT compile or install anything
cd ${CMESS}
cd ..
mkdir build
cd build
cmake3 ${CMESS} -DPYTHON=ON


# python examples are configured, we can test our python wheels now
make pymess-2.7-run-examples
make pymess-3.4-run-examples

