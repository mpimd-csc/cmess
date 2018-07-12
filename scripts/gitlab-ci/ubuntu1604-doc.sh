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

mkdir ${DEPLOY}     || echo 0



apt-get update --yes || true
apt-get install --yes cmake libopenblas-dev libsuitesparse-dev zlib1g-dev libbz2-dev liblzma-dev
apt-get install --yes libarpack2-dev libmatio-dev libsuperlu-dev
apt-get install --yes python-dev python3-dev python-scipy python3-scipy python-setuptools python3-setuptools
apt-get install --yes python-sphinx python-matplotlib zip
easy_install pip
pip install numpydoc
pip install sphinxcontrib-bibtex
cd ..
mkdir -p build
cd build
BUILDDIR=$(pwd)
cmake ${CMESS} -DDEBUG=OFF -DARPACK=ON -DMATIO=ON -DPYTHON=ON -DOPENMP=OFF -DDOC=ON -DDOCPDF=ON
make

PYTHON_VERSION=$(python -c "print '%d.%d'%( __import__('sys').version_info[:2])")
PYTHON3_VERSION=$(python3 -c "print('%d.%d'%( __import__('sys').version_info[:2]))")
make pymess-${PYTHON_VERSION}-build
make pymess-${PYTHON_VERSION}-install
make pymess-${PYTHON3_VERSION}-build
make pymess-${PYTHON3_VERSION}-install

make doc
make pymess-python2.7-doc-html
make pymess-python2.7-doc-pdf

# copy data to deploy directory
cp doc/latex/refman.pdf                     ${DEPLOY}/cmess_doc.pdf
cp python/doc/python2.7/latex/pymess.pdf    ${DEPLOY}/pymess_doc.pdf
zip -r ${DEPLOY}/cmess_doc_html.zip         doc/html/
zip -r ${DEPLOY}/pymess_doc_html.zip        python/doc/python2.7/html/

