#!/bin/sh

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

# call example ./matlab_run.sh matlab2013b

MNAME="$1"
MROOT=$(${MNAME} -n | grep -P "MATLAB\s+=" | awk '{print $5}')

rm -rf mess_build_${MNAME}
mkdir mess_build_${MNAME}
cd mess_build_${MNAME}
#CC=gcc-4.4 CXX=g++-4.4 FC=gfortran-4.4 cmake ../cmess -DDEBUG=OFF -DSUPERLU=OFF -DMATLAB=ON -DMATLAB_ROOT=${MROOT}
cmake ../cmess -DDEBUG=ON -DSUPERLU=OFF -DMATLAB=ON  -DARPACK=ON -DMATLAB_ROOT=${MROOT}
make -j4
${MROOT}/bin/matlab -nodesktop -nojvm -debug  -r "rehash toolboxcache;addpath('$HOME'); run matlab/mess_path.m; ret=run_tests(); exit(ret)"
${MROOT}/bin/matlab -nodesktop -nojvm -debug  -r "rehash toolboxcache;addpath('$HOME'); run matlab/mess_path.m; ret=run_examples(); exit(ret)"

