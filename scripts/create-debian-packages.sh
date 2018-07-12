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

####################################################################
# create-debian-packages.sh generates the debian packages
#
#   - libmess
#   - libmess-doc
#   - libmess-dev
#
# For the distros:
#   - xenial
#   - zesty
#   - jessie
#   - stretch
#
# Execute this script from the root directory of C-M.E.S.S.
#
# Warning it may take several hours until this script has finished.
#
# Unfortunately the pbuilder-dist tool seems to have some
# bugs. The resulting *.deb packages are written to
# to ${HOME}/pbuilder/*.
#
# You have to be able to use sudo because pbuild-dist uses chroot.
#
###################################################################

### VARIABLES FOR CONTROLLING THE SCRIPT    ###

# output directory of packages
PBUILDDIR=$(pwd)/packages
PBUILDERRC=$HOME/.pbuilderrc

# variables for pbuilderrc
APTCACHE=/var/cache/pbuilder/aptcache
AUTOCLEANAPTCACHE=no
http_proxy=http://gitlabci.csc.mpi-magdeburg.mpg.de:3142

#DISTS=(xenial zesty jessie stretch)
DISTS=(xenial zesty jessie stretch)

###############################################

echo "#####################################################################"
echo "# Warning it may take several hours until this script has finished."
echo "# This script produces minimal output."
echo "# Use tail -f on the logfiles created in ${PBUILDDIR}"
echo "# to follow the progress."
echo "# Unfortunately the pbuilder-dist tool seems to have some"
echo "# bugs. The resulting *.deb packages are written to"
echo "# to ${HOME}/pbuilder/* (at least on my machine)."
echo "# You have to be able to use sudo because pbuild-dist uses chroot."
echo "#####################################################################"
sleep 10

###############################################

# make varibles available
export http_proxy=http://gitlabci.csc.mpi-magdeburg.mpg.de:3142

# write to pbuilderrc config
echo "APTCACHE=${APTCACHE}" > ${PBUILDERRC}
echo "AUTOCLEANAPTCACHE=${AUTOCLEANAPTCACHE}" >> ${PBUILDERRC}

# create folder
for DIST in ${DISTS[@]};
    do
    mkdir -p ${PBUILDDIR}/${DIST}
done



# stop on error
set -e -x

# current directory is cmess root directory
CMESS=$(pwd)

# install packages for debian packaging
sudo apt-get update     --yes   > /dev/null
sudo apt-get install    --yes  pbuilder ubuntu-dev-tools debootstrap devscripts > /dev/null

# create environment for pbuilder
for DIST in ${DISTS[@]};
        do
        echo "Current time : $(date +'%T')"
        echo "Dist : ${DIST}"
        pbuilder-dist ${DIST} create > /dev/null
        pbuilder-dist ${DIST} update > /dev/null
done


# generate debian files for pbuilder
echo "Current time : $(date +'%T')"
dpkg-buildpackage  -us -uc > ${PBUILDDIR}/dpkg-buildpackage.log


# create
for DIST in ${DISTS[@]};
        do
        echo "Current time : $(date +'%T')"
        pbuilder-dist ${DIST} ../*.dsc --http-proxy ${http_proxy} --buildresult ${PBUILDDIR}/${DIST} >  ${PBUILDDIR}/${DIST}/pbuilder_${DIST}.log
done





