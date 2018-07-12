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



## This file should be placed in the root directory of your project.
## Then modify the CMakeLists.txt file in the root directory of your
## project to incorporate the testing dashboard.
## # The following are required to uses Dart and the Cdash dashboard
##   ENABLE_TESTING()
##   INCLUDE(CTest)
SET(CTEST_PROJECT_NAME "CMESS")
SET(CTEST_NIGHTLY_START_TIME "01:00:00 UTC")

SET(CTEST_DROP_METHOD "http")
SET(CTEST_DROP_SITE "johannes")
SET(CTEST_DROP_LOCATION "/cdash/submit.php?project=CMESS")
SET(CTEST_DROP_SITE_CDASH TRUE)

#increase timeout for valgrind test to 1500s
SET(DART_TESTING_TIMEOUT "1500")
SET(CTEST_TEST_TIMEOUT "1500")


#valgrind
FIND_PROGRAM(VALGRIND valgrind)
IF(VALGRIND)
    MESSAGE(STATUS "SET valgrind:" ${VALGRIND})

    #increase timeout for valgrind test to 1500s
    SET(DART_TESTING_TIMEOUT "1500")
    SET(CTEST_TEST_TIMEOUT "1500")

    #SET(CTEST_MEMORYCHECK_COMMAND ${VALGRIND})
    SET(MEMORYCHECK_COMMAND ${VALGRIND})

    #SET(CTEST_MEMORYCHECK_COMMAND_OPTIONS "--trace-children=yes --leak-check=full --track-origins=yes")
    SET(MEMORYCHECK_COMMAND_OPTIONS "--trace-children=yes --leak-check=full --track-origins=yes --show-reachable=yes --malloc-fill=0x80 --free-fill=0x7f")
    #SET(CTEST_MEMORYCHECK_SUPPRESSIONS_FILE "${CMAKE_SOURCE_DIR}/cmess_valgrind.supp")
    SET(MEMORYCHECK_SUPPRESSIONS_FILE "${CMAKE_SOURCE_DIR}/cmess_valgrind.supp")

ENDIF()

