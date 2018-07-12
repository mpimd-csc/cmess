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
NUMPY_PY36=http://pc1184.mpi-magdeburg.mpg.de:8080/soft/numpy-1.14.5-cp36-cp36m-manylinux1_x86_64.whl

SCIPY_PY27=http://pc1184.mpi-magdeburg.mpg.de:8080/soft/scipy-0.16.0-cp27-cp27mu-manylinux1_x86_64.whl
SCIPY_PY34=http://pc1184.mpi-magdeburg.mpg.de:8080/soft/scipy-0.16.0-cp34-cp34m-manylinux1_x86_64.whl
SCIPY_PY35=http://pc1184.mpi-magdeburg.mpg.de:8080/soft/scipy-0.16.0-cp35-cp35m-manylinux1_x86_64.whl
SCIPY_PY36=http://pc1184.mpi-magdeburg.mpg.de:8080/soft/scipy-1.1.0-cp36-cp36m-manylinux1_x86_64.whl


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
yum install -y yum-utils
yum install -y centos-release-SCL
yum -y update


# install python2.7, python3.4, python3.5
yum install -y python27 rh-python34 rh-python35 rh-python36


# make python27, python34, python35, python36 available
export PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/rh/devtoolset-2/root/usr/bin
export PATH=/opt/rh/python27/root/usr/bin:$PATH
export PATH=/opt/rh/rh-python34/root/usr/bin:$PATH
export PATH=/opt/rh/rh-python35/root/usr/bin:$PATH
export PATH=/opt/rh/rh-python36/root/usr/bin:$PATH

export LD_LIBRARY_PATH=/usr/local/lib64:/usr/local/lib:/opt/rh/devtoolset-2/root/usr/lib64:/opt/rh/devtoolset-2/root/usr/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/opt/rh/python27/root/usr/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/opt/rh/rh-python34/root/usr/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/opt/rh/rh-python35/root/usr/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/opt/rh/rh-python36/root/usr/lib64:$LD_LIBRARY_PATH


# update pip, install numpy scipy
pip2.7 install --upgrade pip
pip2.7 install ${NUMPY_PY27} ${SCIPY_PY27} wheel
pip3.4 install --upgrade pip
pip3.4 install ${NUMPY_PY34} ${SCIPY_PY34} wheel
pip3.5 install --upgrade pip
pip3.5 install ${NUMPY_PY35} ${SCIPY_PY35} wheel
pip3.6 install --upgrade pip
pip3.6 install ${NUMPY_PY36} ${SCIPY_PY36} wheel



# now install the deployed PYMESS
cd ${DEPLOY}
pip2.7 install *cp27*.whl
pip3.4 install *cp34*.whl
pip3.5 install *cp35*.whl
pip3.6 install *cp36*.whl

# install packages such that we can configure a cmess build
yum install -y cmake3 bzip2-devel zlib-devel xz-devel suitesparse-devel lapack-devel blas-devel


# configure C-M.E.S.S. build and NOT compile or install anything
cd ${CMESS}
cd ..
mkdir build
cd build
cmake3 ${CMESS} -DPYTHON=ON


# python examples are configured, we can test our python wheels now
make pymess-2.7-run-examples
make pymess-3.4-run-examples
make pymess-3.5-run-examples
make pymess-3.6-run-examples


