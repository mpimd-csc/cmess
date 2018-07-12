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

PATCHELF=http://pc1184.mpi-magdeburg.mpg.de:8080/soft/patchelf/patchelf-0.9.tar.gz
OPENBLAS=http://pc1184.mpi-magdeburg.mpg.de:8080/soft/OpenBLAS/v0.2.20.tar.gz
SUITESPARSE=http://pc1184.mpi-magdeburg.mpg.de:8080/soft/SuiteSparse-4.5.6.tar.gz

#additional directories for build install and deploy
DEPBUILD=/depbuild
DEPINSTALL=/depinstall
DEPLOG=/deplog

#check if MESS_ARTIFACTS variable is defined from .gitlab-ci.yml
if [ -z "${MESS_ARTIFACTS}" ]
    then
    echo  "MESS_ARTIFACTS variable is not defined"
    DEPLOY=${CMESS}/mess_artifacts
    else
    echo  "MESS_ARTIFACTS variable is defined"
    DEPLOY=${CMESS}/${MESS_ARTIFACTS}
fi

#check if MESS_VERSION variable is defined from .gitlab-ci.yml
if [ -z "${MESS_VERSION}" ]
    then
    echo  "MESS_VERSION variable is not defined"
    MESS_VERSION="1.0.0"
fi


# strip the binaries
STRIP=false

# allowed library dependencies due to https://www.python.org/dev/peps/pep-0513/
# this string is usefull with grep -vE "${PEP513_GREPPATTERN}"
# note that we have to escape all dots
PEP513_GREPPATTERN="libpanelw\.so\.5"
PEP513_GREPPATTERN="${PEP513_GREPPATTERN}|libncursesw\.so\.5"
PEP513_GREPPATTERN="${PEP513_GREPPATTERN}|libgcc_s\.so\.1"
PEP513_GREPPATTERN="${PEP513_GREPPATTERN}|libstdc++\.so\.6"
PEP513_GREPPATTERN="${PEP513_GREPPATTERN}|libm\.so\.6"
PEP513_GREPPATTERN="${PEP513_GREPPATTERN}|libdl\.so\.2"
PEP513_GREPPATTERN="${PEP513_GREPPATTERN}|librt\.so\.1"
PEP513_GREPPATTERN="${PEP513_GREPPATTERN}|libcrypt\.so\.1"
PEP513_GREPPATTERN="${PEP513_GREPPATTERN}|libc\.so\.6"
PEP513_GREPPATTERN="${PEP513_GREPPATTERN}|libnsl\.so\.1"
PEP513_GREPPATTERN="${PEP513_GREPPATTERN}|libutil\.so\.1"
PEP513_GREPPATTERN="${PEP513_GREPPATTERN}|libpthread\.so\.0"
PEP513_GREPPATTERN="${PEP513_GREPPATTERN}|libresolv\.so\.2"
PEP513_GREPPATTERN="${PEP513_GREPPATTERN}|libX11\.so\.6"
PEP513_GREPPATTERN="${PEP513_GREPPATTERN}|libXext\.so\.6"
PEP513_GREPPATTERN="${PEP513_GREPPATTERN}|libXrender\.so\.1"
PEP513_GREPPATTERN="${PEP513_GREPPATTERN}|libICE\.so\.6"
PEP513_GREPPATTERN="${PEP513_GREPPATTERN}|libSM\.so\.6"
PEP513_GREPPATTERN="${PEP513_GREPPATTERN}|libGL\.so\.1"
PEP513_GREPPATTERN="${PEP513_GREPPATTERN}|libgobject-2\.0\.so\.0"
PEP513_GREPPATTERN="${PEP513_GREPPATTERN}|libgthread-2\.0\.so\.0"
PEP513_GREPPATTERN="${PEP513_GREPPATTERN}|libglib-2\.0\.so\.0"


# create root dirs for build, logging of build, install and deploy
mkdir ${DEPBUILD}
mkdir ${DEPINSTALL}
mkdir ${DEPLOY} || echo 0
mkdir ${DEPLOG}

# update yum
yum clean all 
yum makecache
yum update -y
yum clean all
yum install -y centos-release-SCL
yum -y update

# install packages for C-M.E.S.S.
yum install -y cmake3 bzip2-devel zlib-devel xz-devel

# install packages for patching
yum install -y unzip wget zip curl

# install python2.7, python3.4, python3.5, python3.6
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


# install patchelf
cd ${DEPBUILD}
wget ${PATCHELF} -O patchelf.tar.gz
mkdir patchelf
tar xf patchelf.tar.gz  -C patchelf --strip-components 1
cd patchelf
./configure >> ${DEPLOG}/patchelf.log   2>&1
make >> ${DEPLOG}/patchelf.log          2>&1
make install >> ${DEPLOG}/patchelf.log  2>&1
cd ..
rm -rf patchelf
rm -rf patchelf.tar.gz


# compile openblas on our own
cd ${DEPBUILD}
wget ${OPENBLAS} -O openblas.tar.gz
mkdir openblas
tar xf openblas.tar.gz  -C openblas --strip-components 1
cd openblas
make DYNAMIC_ARCH=1 NO_AVX=1 NO_AVX2=1  >> ${DEPLOG}/openblas.log          2>&1
make PREFIX=${DEPINSTALL} install       >> ${DEPLOG}/openblas.log          2>&1
cd ..
rm -rf openblas
rm -rf openblas.tar.gz


# compile suitesparse on our own
cd ${DEPBUILD}
wget ${SUITESPARSE} -O suitesparse.tar.gz
mkdir suitesparse
tar xf suitesparse.tar.gz  -C suitesparse --strip-components 1
cd suitesparse

# modify SuiteSparse_config.mk
SUITESPARSECONFIGMK=SuiteSparse_config/SuiteSparse_config.mk
sed -i -e 's/CFOPENMP ?=.*/CFOPENMP ?=/g'                                                                                       ${SUITESPARSECONFIGMK}
sed -i -e 's/LDFLAGS +=.*/LDFLAGS += -L$(INSTALL_LIB) -L\'"$DEPINSTALL"'\/lib/g'                                                ${SUITESPARSECONFIGMK}
sed -i -e 's/LAPACK ?=.*/LAPACK ?= -lopenblas/g'                                                                                ${SUITESPARSECONFIGMK}
sed -i -e 's/CF ?=.*/CF ?=$(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) $(OPTIMIZATION) -fexceptions -fPIC -L\'"$DEPINSTALL"'\/lib/g'    ${SUITESPARSECONFIGMK}

export LD_LIBRARY_PATH=/depinstall/lib:$LD_LIBRARY_PATH
make                                    >> ${DEPLOG}/suitesparse.log          2>&1
make library                            >> ${DEPLOG}/suitesparse.log          2>&1
make install INSTALL=${DEPINSTALL}      >> ${DEPLOG}/suitesparse.log          2>&1
cd ..
rm -rf suitesparse
rm -rf suitesparse.tar.gz


# compile arpack on our own
cd ${DEPBUILD}
git clone https://github.com/opencollab/arpack-ng.git
cd arpack-ng
git checkout 3.5.0
mkdir build
cd build
cmake ../ -DBUILD_SHARED_LIBS=ON -DCMAKE_INSTALL_PREFIX=${DEPINSTALL} -DBLAS_LIBRARIES=${DEPINSTALL}/lib/libopenblas.so -DLAPACK_LIBRARIES=${DEPINSTALL}/lib/libopenblas.so -DEXAMPLES=ON
make            >> ${DEPLOG}/arpack.log          2>&1
make install    >> ${DEPLOG}/arpack.log          2>&1
cd ../..
rm -rf arpack-ng


# compile superlu on our own
cd ${DEPBUILD}
git clone https://github.com/xiaoyeli/superlu.git
cd superlu
git checkout v5.2.0
mkdir build
cd build
cmake ../ -Denable_blaslib=OFF -DBLAS_LIBRARIES=/depinstall/lib/libopenblas.so -DCMAKE_INSTALL_PREFIX=/depinstall -DBUILD_SHARED_LIBS=ON
make                >> ${DEPLOG}/superlu.log          2>&1
make install        >> ${DEPLOG}/superlu.log          2>&1


# make a C-M.E.S.S. build
cd /
mkdir -p build
cd build
cmake3 ${CMESS} -DDEBUG=OFF -DARPACK=${DEPINSTALL} -DPYTHON=ON -DOPENMP=OFF -DSUPERLU=ON -DBLAS=${DEPINSTALL}/lib/libopenblas.so \
                -DLAPACK=${DEPINSTALL}/lib/libopenblas.so -DSUITESPARSE=${DEPINSTALL} -DSUPERLU_ROOT=${DEPINSTALL}
make

make pymess-2.7-build
make pymess-3.4-build
make pymess-3.5-build
make pymess-3.6-build

################################# DESCRIPTION FOR WHEEL CREATION ################################################################
#                                                                                                                               #
# Change to Python Build Directory and create bdist_wheel in the directory dist.                                                #
# Unzip the wheel and remove the wheel.                                                                                         #
# Change to unzipped pymess directory.                                                                                          #
#                                                                                                                               #
# Get Dependencies of _c_interface*.so*, copy them to ./libs directory and set the rpath for each lib in ./libs:                #
#   1.  grep for =>                                                                                                             #
#   2.  remove content in brackets version ect                                                                                  #
#   3.  remove leading whitespaces                                                                                              #
#   4.  remove double whitespaces                                                                                               #
#   5.  remove trailing whitespaces                                                                                             #
#   6.  remove all shown dependencies with no target (String ends with "=>")                                                    #
#   7.  remove double entries with sort | uniq                                                                                  #
#   8.  remove allowed library dependencies due to https://www.python.org/dev/peps/pep-0513/                                    #
#   9.  use awk to split between "=>" ( <name of library> => <absolute path of library> ),                                      #
#       copy then the library to ./libs (copy is used with overwrite)                                                           #
#   10. set RPATH to $ORIGIN for each library in ./libs,                                                                        #
#       (the RPATH of _c_interface*.so points already to $ORIGIN/.libs, see setup.py)                                           #
#   11. move back and zip together                                                                                              #
#   12. copy everything to ${DEPLOY}                                                                                            #
#                                                                                                                               #
# Move upwards and zip the modified archives to a wheel and copy everything to deploy folder.                                   #
#                                                                                                                               #
#################################################################################################################################


#------------------------------ CREATE A WHEEL FOR PYTHON2.7  ---------------------------------------------#
cd /build/python/python_2.7
python2.7 setup.py bdist_wheel
cd dist
unzip *.whl
rm -rf *.whl
cd pymess

ldd -v _c_interface*.so*    | grep "=>" | sed -e 's/([^()]*)//g' | sed -e 's/^[ \t]*//' | sed 's/  */ /g' | sed 's/\s*$//g'     \
                            | grep -vE "=>$" | sort  | uniq | grep -vE "${PEP513_GREPPATTERN}"                                  \
                            | awk -F '=>' '{print "PERFORM PATCH COMMAND:  cp "$2" ./libs/"$1""; system("cp "$2" .libs/"$1""); }'

for lib in $(find .libs -name "*.so*");
    do
    librpath=$(patchelf --print-rpath ${lib});
    echo "Set RPATH of ${lib} from ${librpath} to '\$ORIGIN'";
    patchelf --set-rpath '$ORIGIN' $lib;
done

if test "$STRIP" = true; then
    echo "STRIP ALL LIBS"
    for lib in $(find .libs -name "*.so*");
        do
        echo "Strip ${lib}";
        strip -s $lib;
    done
fi

cd ..
zip -r pymess-${MESS_VERSION}-cp27-cp27mu-manylinux1_x86_64.whl  .
cp  pymess-${MESS_VERSION}-cp27-cp27mu-manylinux1_x86_64.whl  ${DEPLOY}
############################################################################################################


#------------------------------ CREATE A WHEEL FOR PYTHON3.4  ---------------------------------------------#
cd /build/python/python_3.4
python3.4 setup.py bdist_wheel
cd dist
unzip *.whl
rm -rf *.whl
cd pymess

ldd -v _c_interface*.so*    | grep "=>" | sed -e 's/([^()]*)//g' | sed -e 's/^[ \t]*//' | sed 's/  */ /g' | sed 's/\s*$//g'     \
                            | grep -vE "=>$" | sort  | uniq | grep -vE "${PEP513_GREPPATTERN}"                                  \
                            | awk -F '=>' '{print "PERFORM PATCH COMMAND:  cp "$2" ./libs/"$1""; system("cp "$2" .libs/"$1""); }'

for lib in $(find .libs -name "*.so*");
    do
    librpath=$(patchelf --print-rpath ${lib});
    echo "Set RPATH of ${lib} from ${librpath} to '\$ORIGIN'";
    patchelf --set-rpath '$ORIGIN' $lib;
done

if test "$STRIP" = true; then
    echo "STRIP ALL LIBS"
    for lib in $(find .libs -name "*.so*");
        do
        echo "Strip ${lib}";
        strip -s $lib;
    done
fi

cd ..
zip -r pymess-${MESS_VERSION}-cp34-cp34m-manylinux1_x86_64.whl  .
cp  pymess-${MESS_VERSION}-cp34-cp34m-manylinux1_x86_64.whl  ${DEPLOY}
############################################################################################################



#------------------------------ CREATE A WHEEL FOR PYTHON3.5  ---------------------------------------------#
cd /build/python/python_3.5
python3.5 setup.py bdist_wheel
cd dist
unzip *.whl
rm -rf *.whl
cd pymess

ldd -v _c_interface*.so*    | grep "=>" | sed -e 's/([^()]*)//g' | sed -e 's/^[ \t]*//' | sed 's/  */ /g' | sed 's/\s*$//g'     \
                            | grep -vE "=>$" | sort  | uniq | grep -vE "${PEP513_GREPPATTERN}"                                  \
                            | awk -F '=>' '{print "PERFORM PATCH COMMAND:  cp "$2" ./libs/"$1""; system("cp "$2" .libs/"$1""); }'

for lib in $(find .libs -name "*.so*");
    do
    librpath=$(patchelf --print-rpath ${lib});
    echo "Set RPATH of ${lib} from ${librpath} to '\$ORIGIN'";
    patchelf --set-rpath '$ORIGIN' $lib;
done

if test "$STRIP" = true; then
    echo "STRIP ALL LIBS"
    for lib in $(find .libs -name "*.so*");
        do
        echo "Strip ${lib}";
        strip -s $lib;
    done
fi

cd ..
zip -r pymess-${MESS_VERSION}-cp35-cp35m-manylinux1_x86_64.whl  .
cp  pymess-${MESS_VERSION}-cp35-cp35m-manylinux1_x86_64.whl  ${DEPLOY}
############################################################################################################



#------------------------------ CREATE A WHEEL FOR PYTHON3.6  ---------------------------------------------#
cd /build/python/python_3.6
python3.6 setup.py bdist_wheel
cd dist
unzip *.whl
rm -rf *.whl
cd pymess

ldd -v _c_interface*.so*    | grep "=>" | sed -e 's/([^()]*)//g' | sed -e 's/^[ \t]*//' | sed 's/  */ /g' | sed 's/\s*$//g'     \
                            | grep -vE "=>$" | sort  | uniq | grep -vE "${PEP513_GREPPATTERN}"                                  \
                            | awk -F '=>' '{print "PERFORM PATCH COMMAND:  cp "$2" ./libs/"$1""; system("cp "$2" .libs/"$1""); }'

for lib in $(find .libs -name "*.so*");
    do
    librpath=$(patchelf --print-rpath ${lib});
    echo "Set RPATH of ${lib} from ${librpath} to '\$ORIGIN'";
    patchelf --set-rpath '$ORIGIN' $lib;
done

if test "$STRIP" = true; then
    echo "STRIP ALL LIBS"
    for lib in $(find .libs -name "*.so*");
        do
        echo "Strip ${lib}";
        strip -s $lib;
    done
fi

cd ..
zip -r pymess-${MESS_VERSION}-cp36-cp36m-manylinux1_x86_64.whl  .
cp  pymess-${MESS_VERSION}-cp36-cp36m-manylinux1_x86_64.whl  ${DEPLOY}
############################################################################################################



# SOME TRIES WHICH NOT REALLY WORKED OR ARE TO COMPLICATED

######################## COPY LIBRARIES ###############################
# 1. We copy the needed libraries files manually to .libs.
# 2. We get the needed libraries from
#    readelf -d libmess_python_2.7.so | grep 'NEEDED'
# 3. The corresponding Path we get from ldd libmess_python_2.7.so
# 4. The Last step is to set the RPATH to $ORIGIN for all libraries
#    and check again the dependencies with `readelf ...`
#
#for LIB in  libamd.so libmetis.so libcamd.so libcolamd.so libcxsparse.so libopenblas.so     \
#            libsuitesparseconfig.so libsuperlu.so libumfpack.so libarpack.so;
#    do
#    find ${DEPINSTALL} -name "${LIB}" -exec cp {} . \;
#done
#cp /usr/lib64/libbz2.so                             .
#cp /usr/lib64/liblzma.so                            .
#cp /usr/lib64/libz.so                               .
#cp /opt/rh/python27/root/usr/lib64/libpython2.7.so  .
#
#
## patch amd
#patchelf --replace-needed 'libsuitesparseconfig.so.4' 'libsuitesparseconfig.so' libamd.so
#
## patch arpack
#patchelf --replace-needed 'libopenblas.so.0' 'libopenblas.so' libarpack.so
#
## patch camd
#patchelf --replace-needed 'libsuitesparseconfig.so.4'       'libsuitesparseconfig.so'    libcamd.so
#
## patch cholmod
#patchelf --replace-needed 'libamd.so.2'                     'libamd.so'                 libcholmod.so
#patchelf --replace-needed 'libcolamd.so.2'                  'libcolamd.so'              libcholmod.so
#patchelf --replace-needed 'libsuitesparseconfig.so.4'       'libsuitesparseconfig.so'   libcholmod.so
#patchelf --replace-needed 'libccolamd.so.2'                 'libccolamd.so'             libcholmod.so
#patchelf --replace-needed 'libcamd.so.2'                    'libcamd.so'                libcholmod.so
#patchelf --replace-needed 'libopenblas.so.0'                'libopenblas.so'            libcholmod.so
#
## patch colamd
#patchelf --replace-needed 'libsuitesparseconfig.so.4'       'libsuitesparseconfig.so'   libcolamd.so
#
## patch cxsparse
#patchelf --replace-needed 'libsuitesparseconfig.so.4'       'libsuitesparseconfig.so'   libcxsparse.so
#
## patch pymess
#patchelf --replace-needed 'libsuitesparseconfig.so.4'       'libsuitesparseconfig.so'   libmess_python_2.7.so
#patchelf --replace-needed 'libopenblas.so.0'                'libopenblas.so'            libmess_python_2.7.so
#patchelf --replace-needed 'libamd.so.2'                     'libamd.so'                 libmess_python_2.7.so
#patchelf --replace-needed 'libcolamd.so.2'                  'libcolamd.so'              libmess_python_2.7.so
#patchelf --replace-needed 'libumfpack.so.5'                 'libumfpack.so'             libmess_python_2.7.so
#patchelf --replace-needed 'libcholmod.so.3'                 'libcholmod.so'             libmess_python_2.7.so
#patchelf --replace-needed 'libcxsparse.so.3'                'libcxsparse.so'            libmess_python_2.7.so
#patchelf --replace-needed 'libsuperlu.so.5'                 'libsuperlu.so'             libmess_python_2.7.so
#patchelf --replace-needed 'libarpack.so.2'                  'libarpack.so'              libmess_python_2.7.so
#patchelf --replace-needed 'liblzma.so.0'                    'liblzma.so'                libmess_python_2.7.so
#patchelf --replace-needed 'libbz2.so.1'                     'libbz2.so'                 libmess_python_2.7.so
#patchelf --replace-needed 'libz.so.1'                       'libz.so'                   libmess_python_2.7.so
#
## patch superlu
#patchelf --replace-needed 'libopenblas.so.0'                'libopenblas.so'            libsuperlu.so
#
## patch umfpack
#patchelf --replace-needed 'libamd.so.2'                     'libamd.so'                 libumfpack.so
#patchelf --replace-needed 'libsuitesparseconfig.so.4'       'libsuitesparseconfig.so'   libumfpack.so
#patchelf --replace-needed 'libcholmod.so.3'                 'libcholmod.so'             libumfpack.so
#
#find . -name "*.so*" -exec patchelf --set-rpath '$ORIGIN' {} \;
## patch _c_interface
#cd ..
#CINTERFACE=$(find . -name "_c_interface*" -exec basename {} \;)
#patchelf --replace-needed 'libpython2.7.so.1.0'             'libpython2.7.so'           ${CINTERFACE}
#
#
##### MOVE UP AND ZIP
#cd ..
#zip -r pymess-0.9.0-cp27-cp27mu-manylinux1_x86_64.whl  .
#cp  pymess-0.9.0-cp27-cp27mu-manylinux1_x86_64.whl  /deploy/
#
################################################################################################################


####################################### GET DEPENDENCIES################################################################
# 1. get dependencies of libpymess*.so with readelf
# 2. take only dependencies marked with NEEDED
# 3. use awk to get the library names in brackets
# 4. remove the libraries which we do not need to ship see also https://www.python.org/dev/peps/pep-0513/
# 5. store the resulting libs in an array
#
# 1. iterate over the library dependencies
# 2. get the absolute path from ldd
# 3. copy the library to that place (pymess/.libs)
#
# this string is usefull with grep -vE "${PEP513_GREPPATTERN}"

#DEPLIBS_=$(readelf -d *.so | grep "NEEDED" | awk -F'[][]' '{print $2}' | grep -vE "${PEP513_GREPPATTERN}")
#for DEPLIB_ in ${DEPLIBS_};
#    do
#    echo ${DEPLIB_}
#
#    done
#######################################################################################################################


