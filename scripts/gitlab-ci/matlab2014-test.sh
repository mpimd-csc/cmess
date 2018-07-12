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
apt-get install --yes gcc-4.9 g++-4.9 gfortran-4.9
cd ..
mkdir -p build
cd build

#get matlab root
MROOT=$(matlab -n | grep -P "MATLAB\s+=" | awk '{print $5}')

#compile mex-mess
CC=gcc-4.9 CXX=g++-4.9 FC=gfortran-4.9 cmake ${CMESS} -DDEBUG=OFF -DSUPERLU=OFF -DARPACK=ON -DOPENMP=OFF -DMATLAB=ON -DMATLAB_ROOT=${MROOT}
make
export LD_LIBRARY_PATH=$(pwd)/lib/:$LD_LIBRARY_PATH

#run mex-mess tests
${MROOT}/bin/matlab -nodesktop -nojvm -debug -r "rehash toolboxcache; addpath('$HOME'); run matlab/mess_path.m; ret=run_tests(); exit(ret)"
${MROOT}/bin/matlab -nodesktop -nojvm -debug -r "rehash toolboxcache; addpath('$HOME'); run matlab/mess_path.m; ret=run_examples(); exit(ret)"
