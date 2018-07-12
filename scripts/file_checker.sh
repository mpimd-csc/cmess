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

# file_checker.sh searches for
#
#   * Windows Carriage Return
#   * Tab Characters
#   * Trailing Whitespaces
#



echo "##################### FILE-CHECKER: Search For Windows Carriage Return #############################"
WCR=$(grep      -o  --exclude=\*.{pyc,png,bib,mod,ico,gitmodules} --exclude-dir=csc-bibfiles --exclude-dir=obj-x86_64-linux-gnu  -rn $'\r' *)
WCR_CNT=$(grep  -o  --exclude=\*.{pyc,png,bib,mod,ico,gitmodules} --exclude-dir=csc-bibfiles --exclude-dir=obj-x86_64-linux-gnu  -rn $'\r' * |  wc -l)

echo "$WCR"
echo "Found: $WCR_CNT"

if [[ ${WCR_CNT} -ge 1 ]]; then exit 1; fi
echo -e "-----------------------------------------------------------------------------------------------------\n\n"

echo "##################### FILE-CHECKER: Search For Tab Character #######################################"
TABS=$(grep     -o --exclude=\*.{pyc,png,bib,mod,ico,gitmodules} --exclude=tags --exclude-dir=csc-bibfiles --exclude-dir=obj-x86_64-linux-gnu  --exclude=\*rules   -rn $'\t' *)
TABS_CNT=$(grep -o --exclude=\*.{pyc,png,bib,mod,ico,gitmodules} --exclude=tags --exclude-dir=csc-bibfiles --exclude-dir=obj-x86_64-linux-gnu  --exclude=\*rules   -rn $'\t' * | wc -l)

echo "$TABS"
echo "Found: $TABS_CNT"

if [[ ${TABS_CNT} -ge 1 ]]; then exit 1; fi
echo -e "-----------------------------------------------------------------------------------------------------\n\n"

echo "##################### FILE-CHECKER: Search For Trailing Whitespaces ################################"
TRAILINGWHITESPACES=$(grep     -o --exclude=\*.{pyc,png,bib,mod,ico,gitmodules} --exclude-dir=csc-bibfiles  --exclude-dir=obj-x86_64-linux-gnu  -rnL '[^[:blank:]]$' *)
TRAILINGWHITESPACES_CNT=$(grep -o --exclude=\*.{pyc,png,bib,mod,ico,gitmodules} --exclude-dir=csc-bibfiles  --exclude-dir=obj-x86_64-linux-gnu  -rnL '[^[:blank:]]$' * | wc -l)

echo "$TRAILINGWHITESPACES"
echo "Found: $TRAILINGWHITESPACES_CNT"

if [[ ${TRAILINGWHITESPACES_CNT} -ge 1 ]]; then exit 1; fi
echo -e "-----------------------------------------------------------------------------------------------------\n\n"

