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


# get sources for cmess in variables CMESS_SRC and CMESS_FSRC
ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/lib/contrib        ${CMAKE_BINARY_DIR}/python/lib/contrib)
ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/lib/direct         ${CMAKE_BINARY_DIR}/python/lib/direct)
ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/lib/dynsys         ${CMAKE_BINARY_DIR}/python/lib/dynsys)
ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/lib/easyfrontend   ${CMAKE_BINARY_DIR}/python/lib/easyfrontend)
ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/lib/eigen          ${CMAKE_BINARY_DIR}/python/lib/eigen)
ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/lib/extblas        ${CMAKE_BINARY_DIR}/python/lib/extblas)
ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/lib/formats        ${CMAKE_BINARY_DIR}/python/lib/formats)
ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/lib/graph          ${CMAKE_BINARY_DIR}/python/lib/graph)
ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/lib/io             ${CMAKE_BINARY_DIR}/python/lib/io)
ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/lib/itsolver       ${CMAKE_BINARY_DIR}/python/lib/itsolver)
ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/lib/lrcf_adi       ${CMAKE_BINARY_DIR}/python/lib/lrcf_adi)
ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/lib/matrix         ${CMAKE_BINARY_DIR}/python/lib/matrix)
ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/lib/misc           ${CMAKE_BINARY_DIR}/python/lib/misc)
ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/lib/plot           ${CMAKE_BINARY_DIR}/python/lib/plot)
ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/lib/precond        ${CMAKE_BINARY_DIR}/python/lib/precond)
ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/lib/reorder        ${CMAKE_BINARY_DIR}/python/lib/reorder)
ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/lib/vector         ${CMAKE_BINARY_DIR}/python/lib/vector)


# set python version to search for
SET(PY_LIST 2.7 3.0 3.1 3.2 3.3 3.4 3.5 3.6)


IF(DEFINED PYTHON_VERSION)
    SET(PYTHON_VERSION ${PY_VERSION})
ENDIF()


FILE(GLOB_RECURSE DEP_C "*.[hc]")
FILE(GLOB_RECURSE DEP_CC "*.cc")
FILE(GLOB_RECURSE DEP_PY "*.py")
SET(DEP ${DEP_C} ${DEP_CC} ${DEP_PY})



FOREACH(PY ${PY_LIST})
    MESSAGE(STATUS "")
    MESSAGE(STATUS "Look for Python ${PY}")
    SET(PythonInterp_FIND_VERSION ${PY})

    #set the desired version to look for
    SET(Python_ADDITIONAL_VERSIONS ${PY})

    FIND_PACKAGE(PythonInterp)

    IF("${PY}" STREQUAL "${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}")
        SET(NUMPY_FOUND FALSE)
        SET(SCIPY_FOUND FALSE)
        FIND_PACKAGE(NumPy)
        FIND_PACKAGE(SciPy)

        IF(NUMPY_FOUND AND SCIPY_FOUND)

            # update vars for configure
            SET(PYMESS_OUTPUT_DIR           "${CMAKE_CURRENT_BINARY_DIR}/python_${PY}")
            SET(PYMESS_PYTHON_LIB           mess_python_${PY})
            SET(PYMESS_LIBRARY_LIST         "'${PYMESS_PYTHON_LIB}'")
            SET(PYMESS_LIBRARY_PATH_LIST    "'${PYMESS_OUTPUT_DIR}/pymess/.libs'")
            #SET(PYMESS_LIBRARY_PATH_LIST    "'${PYMESS_OUTPUT_DIR}/pymess/'")

            # add a special cmess library for this python version
            ADD_LIBRARY(${PYMESS_PYTHON_LIB}  ${CMESS_SRC} ${CMESS_FSRC})
            TARGET_INCLUDE_DIRECTORIES(${PYMESS_PYTHON_LIB} PUBLIC ${MESS_INCLUDE_DIR} ${CMAKE_BINARY_DIR} ${PYTHON_INCLUDE_DIR})
            TARGET_LINK_LIBRARIES(${PYMESS_PYTHON_LIB} ${MESS_LIBS})
            SET_TARGET_PROPERTIES(${PYMESS_PYTHON_LIB} PROPERTIES COMPILE_FLAGS "-DPYMESS")
            SET_TARGET_PROPERTIES(${PYMESS_PYTHON_LIB}
                PROPERTIES
                ARCHIVE_OUTPUT_DIRECTORY "${PYMESS_OUTPUT_DIR}/pymess/.libs"
                LIBRARY_OUTPUT_DIRECTORY "${PYMESS_OUTPUT_DIR}/pymess/.libs"
                #ARCHIVE_OUTPUT_DIRECTORY "${PYMESS_OUTPUT_DIR}/pymess/"
                #LIBRARY_OUTPUT_DIRECTORY "${PYMESS_OUTPUT_DIR}/pymess/"
            )


            # configure setup.py
            SET(PYMESS_INCLUDE_LIST
                "'${CMAKE_SOURCE_DIR}/include',\
                '${CMAKE_BINARY_DIR}/include',\
                '${NUMPY_INCLUDE_DIRS}',\
                '${SCIPY_INCLUDE_DIRS}',\
                '${CMAKE_SOURCE_DIR}/libcscutils/include'")
            SET(SETUP_PY "${PYMESS_OUTPUT_DIR}/setup.py")
            CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in" ${SETUP_PY})


            # copy pymess python sources, because setuptools only allows relative paths
            FILE(GLOB_RECURSE PYMESS_PY_SRCS "${CMAKE_CURRENT_SOURCE_DIR}/pymess/*.py")
            FOREACH(PYMESS_PY_SRC_IN ${PYMESS_PY_SRCS})
                STRING(REPLACE ${CMAKE_CURRENT_SOURCE_DIR} ${PYMESS_OUTPUT_DIR} PYMESS_PY_SRC_OUT ${PYMESS_PY_SRC_IN})
                CONFIGURE_FILE(${PYMESS_PY_SRC_IN} ${PYMESS_PY_SRC_OUT} COPYONLY)
            ENDFOREACH()


            # copy COPYING FILE to build sources
            CONFIGURE_FILE("${CMAKE_SOURCE_DIR}/COPYING" "${PYMESS_OUTPUT_DIR}/COPYING" COPYONLY)


            # copy pymess C sources, because setuptools only allows relative paths
            FILE(GLOB_RECURSE PYMESS_C_SRCS "${CMAKE_CURRENT_SOURCE_DIR}/c_interface/*.c" "${CMAKE_CURRENT_SOURCE_DIR}/c_interface/*.h")
            FOREACH(PYMESS_C_SRC_IN ${PYMESS_C_SRCS})
                STRING(REPLACE ${CMAKE_CURRENT_SOURCE_DIR} "${PYMESS_OUTPUT_DIR}" PYMESS_C_SRC_OUT ${PYMESS_C_SRC_IN})
                CONFIGURE_FILE(${PYMESS_C_SRC_IN} ${PYMESS_C_SRC_OUT} COPYONLY)
            ENDFOREACH()


            # configure / copy examples
            SET(RUN_PY_IN   "${CMAKE_CURRENT_SOURCE_DIR}/example/run.sh.in")
            SET(RUN_PY_OUT  "${PYMESS_OUTPUT_DIR}/pymess/example/run.sh")
            CONFIGURE_FILE(${RUN_PY_IN}  ${RUN_PY_OUT} @ONLY)
            FILE(GLOB_RECURSE PYMESS_EXAMPLES "${CMAKE_SOURCE_DIR}/python/example/*.py")
            FOREACH(PYMESS_EXAMPLE_IN ${PYMESS_EXAMPLES})
                STRING(REPLACE ${CMAKE_CURRENT_SOURCE_DIR} "${PYMESS_OUTPUT_DIR}/pymess" PYMESS_EXAMPLE_OUT ${PYMESS_EXAMPLE_IN})
                CONFIGURE_FILE(${PYMESS_EXAMPLE_IN} ${PYMESS_EXAMPLE_OUT} @ONLY)
            ENDFOREACH()


            # configure test suite
            FILE(GLOB_RECURSE UNITTESTS ${CMAKE_CURRENT_SOURCE_DIR}/unittests/*.py)
            FOREACH(UNITTEST_IN ${UNITTESTS})
                STRING(REPLACE ${CMAKE_CURRENT_SOURCE_DIR} "${PYMESS_OUTPUT_DIR}/pymess" UNITTEST_OUT ${UNITTEST_IN})
                CONFIGURE_FILE(${UNITTEST_IN} ${UNITTEST_OUT} @ONLY)
            ENDFOREACH()

            # configure coverage file and add target
            MESSAGE("")
            MESSAGE_COLOR(STATUS COLOR_BOLD_GREEN "Search Python coverage-${PY} programm.")
            FIND_PROGRAM(PYTHON_COVERAGE NAMES coverage-${PY})

            IF(PYTHON_COVERAGE)
                MESSAGE(STATUS "Found coverage-${PY} programm: ${PYTHON_COVERAGE}")
                MESSAGE_COLOR(STATUS COLOR_BOLD_CYAN "Activate New Target: make pymess-${PY}-coverage")

                CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/python/coverage.rc.in  ${PYMESS_OUTPUT_DIR}/coverage.rc @ONLY)
                ADD_CUSTOM_TARGET("pymess-${PY}-coverage"
                    COMMAND ${PYTHON_COVERAGE} erase
                    COMMAND ${PYTHON_COVERAGE} run      --rcfile=${PYMESS_OUTPUT_DIR}/coverage.rc  ${PYMESS_OUTPUT_DIR}/pymess/unittests/run.py
                    COMMAND ${PYTHON_COVERAGE} report   --rcfile=${PYMESS_OUTPUT_DIR}/coverage.rc
                    COMMAND ${PYTHON_COVERAGE} html     --rcfile=${PYMESS_OUTPUT_DIR}/coverage.rc
                )
            ENDIF()
            UNSET(PYTHON_COVERAGE CACHE)


            # add custom target for setup.py build
            SET(OUTPUT "${PYMESS_OUTPUT_DIR}/build/timestamp-${PY}")
            ADD_CUSTOM_COMMAND(OUTPUT ${OUTPUT}
                COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY} build --force
                COMMAND ${CMAKE_COMMAND} -E touch ${OUTPUT} DEPENDS ${DEP}
                WORKING_DIRECTORY ${PYMESS_OUTPUT_DIR}
                )
            MESSAGE_COLOR(STATUS COLOR_BOLD_CYAN "Activate New Target: make pymess-${PY}-build")
            ADD_CUSTOM_TARGET(pymess-${PY}-build ALL DEPENDS ${OUTPUT})
            ADD_DEPENDENCIES(pymess-${PY}-build ${PYMESS_PYTHON_LIB})


            # add custom target for setup.py install
            MESSAGE_COLOR(STATUS COLOR_BOLD_CYAN "Activate New Target: make pymess-${PY}-install")
            ADD_CUSTOM_TARGET(pymess-${PY}-install
                COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY} install
                WORKING_DIRECTORY ${PYMESS_OUTPUT_DIR}
            )

            # add custom target for run examples
            MESSAGE_COLOR(STATUS COLOR_BOLD_CYAN "Activate New Target: make pymess-${PY}-run-examples")
            ADD_CUSTOM_TARGET(pymess-${PY}-run-examples
                    COMMAND ./run.sh
                    WORKING_DIRECTORY ${PYMESS_OUTPUT_DIR}/pymess/example
                    )

            # add custom target for run tests
            MESSAGE_COLOR(STATUS COLOR_BOLD_CYAN "Activate New Target: make pymess-${PY}-run-tests")
            ADD_CUSTOM_TARGET(pymess-${PY}-run-tests
                    COMMAND ${PYTHON_EXECUTABLE} run.py
                    WORKING_DIRECTORY ${PYMESS_OUTPUT_DIR}/pymess/unittests
                    )


            LIST(APPEND _PYCMESS ${PY})

            # status message that pymess is built
            MESSAGE_COLOR(STATUS COLOR_BOLD_GREEN "Py-M.E.S.S. is built for Python ${PY}: ${PYTHON_EXECUTABLE}")

            # save all valid python executables
            ADD(PYTHON_EXECUTABLES ${PYTHON_EXECUTABLE})

        ELSE()
            MESSAGE(STATUS "Py-M.E.S.S. not build for Python ${PY} because atleast one of NumPy or Scipy is missing ")
        ENDIF()
    ENDIF()
ENDFOREACH()


# MAKE LIST OF AVAILABLE PYTHON VERSIONS AVAILABLE FOR PARENT
SET(PYCMESS ${_PYCMESS} PARENT_SCOPE)




# search tools for pymess documentation target
MESSAGE("")
MESSAGE_COLOR(STATUS COLOR_BOLD_GREEN "Search Python Sphinx Tools for Documentation")
FIND_PACKAGE(SPHINXBUILD)
IF(SPHINXBUILD_FOUND)

    # configure doc conf file
    CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/python/doc/conf.py.in  ${CMAKE_BINARY_DIR}/python/conf.py @ONLY)

    FOREACH(PYTHON_EXECUTABLE ${PYTHON_EXECUTABLES})

        # unset old variables
        UNSET(PY_SPHINXCONTRIB.BIBTEX_FOUND CACHE)
        UNSET(PY_SPHINXCONTRIB.BIBTEX CACHE)
        UNSET(PY_NUMPYDOC_FOUND CACHE)
        UNSET(PY_NUMPYDOC CACHE)
        UNSET(PY_MATPLOTLIB_FOUND CACHE)
        UNSET(PY_MATPLOTLIB CACHE)

        # find a python module, needs PYTHON_EXECUTABLE
        MESSAGE(STATUS "Search Python Sphinx Libraries for: ${PYTHON_EXECUTABLE}")
        FIND_PYTHON_MODULE(${PYTHON_EXECUTABLE} sphinxcontrib.bibtex)
        FIND_PYTHON_MODULE(${PYTHON_EXECUTABLE} numpydoc)
        FIND_PYTHON_MODULE(${PYTHON_EXECUTABLE} matplotlib)

        IF(PY_NUMPYDOC AND PY_SPHINXCONTRIB.BIBTEX AND PY_MATPLOTLIB)
            GET_FILENAME_COMPONENT(PYTHON_EXECUTABLE_NAME ${PYTHON_EXECUTABLE} NAME)

            # setup html output
            SET(TARGET_NAME_HTML "pymess-${PYTHON_EXECUTABLE_NAME}-doc-html")
            SET(PYMESS_DOC_HTML_OUTPUT      "${CMAKE_BINARY_DIR}/python/doc/${PYTHON_EXECUTABLE_NAME}/html")
            SET(PYMESS_CONF_PY_DIR          "${CMAKE_BINARY_DIR}/python/")
            SET(PYMESS_DOC_DIR              "${CMAKE_SOURCE_DIR}/python/doc")
            MESSAGE_COLOR(STATUS COLOR_BOLD_CYAN "Activate New Target: make ${TARGET_NAME_HTML}")
            ADD_CUSTOM_TARGET("${TARGET_NAME_HTML}"
                COMMAND ${SPHINXBUILD} -n -T -v -b doctest  -c ${PYMESS_CONF_PY_DIR} ${PYMESS_DOC_DIR} ${PYMESS_DOC_HTML_OUTPUT}
                COMMAND ${SPHINXBUILD} -n -T -v -b html     -c ${PYMESS_CONF_PY_DIR} ${PYMESS_DOC_DIR} ${PYMESS_DOC_HTML_OUTPUT}
                COMMAND ${CMAKE_COMMAND} -E cmake_echo_color    --cyan "Py-M.E.S.S. html-documentation available: ${PYMESS_DOC_HTML_OUTPUT}"
                COMMENT "Generate Py-M.E.S.S. (html) documentation using ${PYTHON_EXECUTABLE_NAME}"
                )

            # setup latex output
            SET(TARGET_NAME_PDF  "pymess-${PYTHON_EXECUTABLE_NAME}-doc-pdf")
            SET(PYMESS_DOC_LATEX_OUTPUT     "${CMAKE_BINARY_DIR}/python/doc/${PYTHON_EXECUTABLE_NAME}/latex")
            MESSAGE_COLOR(STATUS COLOR_BOLD_CYAN "Activate New Target: make ${TARGET_NAME_PDF}")
            ADD_CUSTOM_TARGET("${TARGET_NAME_PDF}"
                COMMAND ${SPHINXBUILD} -n -T -v -b doctest  -c ${PYMESS_CONF_PY_DIR} ${PYMESS_DOC_DIR} ${PYMESS_DOC_LATEX_OUTPUT}
                COMMAND ${SPHINXBUILD} -n -T -v -b latex    -c ${PYMESS_CONF_PY_DIR} ${PYMESS_DOC_DIR} ${PYMESS_DOC_LATEX_OUTPUT}
                COMMAND make -C ${PYMESS_DOC_LATEX_OUTPUT}  LATEXMKOPTS="-silent"
                COMMAND ${CMAKE_COMMAND} -E cmake_echo_color    --cyan "Py-M.E.S.S. pdf-documentation available: ${PYMESS_DOC_LATEX_OUTPUT}"
                COMMENT "Generate Py-M.E.S.S. (pdf) documentation using ${PYTHON_EXECUTABLE_NAME}"
                )
        ENDIF()
        MESSAGE("")
    ENDFOREACH()
ENDIF()


# search for python code checker: pylint
MESSAGE("")
MESSAGE_COLOR(STATUS COLOR_BOLD_GREEN "Search Python Code Checker: pylint")
FIND_PACKAGE(PYLINT)
IF(PYLINT_FOUND)
    MESSAGE(STATUS "Version: ${PYLINT_VERSION_STRING}")
    IF(${PYLINT_VERSION_STRING} VERSION_LESS "1.6.5")
        MESSAGE(STATUS "Cannot use ${PYLINT}, need at least Version 1.6.5")
    ELSE()
        MESSAGE_COLOR(STATUS COLOR_BOLD_CYAN "Activate New Target: make pymess-lint")
        ADD_CUSTOM_TARGET("pymess-lint"
            COMMAND ${PYLINT} --rcfile=${CMAKE_SOURCE_DIR}/python/pylint.rc --reports=y ${DEP_PY}
            )
    MESSAGE("")
    ENDIF()
ENDIF()




