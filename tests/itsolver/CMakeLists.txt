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

SET(SRC    check_itsolver
    check_precond
    )


FOREACH(srcfile ${SRC})
    ADD_EXECUTABLE(${srcfile} ${srcfile}.c)
    TARGET_LINK_LIBRARIES(${srcfile} ${MESS_LIBS})
ENDFOREACH()


# Add tests
ADD_TEST(check_itsolver_cg ${EXECUTABLE_OUTPUT_PATH}/check_itsolver CG "${CMAKE_SOURCE_DIR}/tests/data/SPD_SYM_625.mtx" 1000 1e-7)
ADD_TEST(check_itsolver_bicgstab ${EXECUTABLE_OUTPUT_PATH}/check_itsolver BICGSTAB "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" 1000 1e-7)
ADD_TEST(check_itsolver_bicgstab_cpx ${EXECUTABLE_OUTPUT_PATH}/check_itsolver BICGSTAB "${CMAKE_SOURCE_DIR}/tests/data/randc50x50.mtx" 1000 1e-7)
ADD_TEST(check_itsolver_gmres ${EXECUTABLE_OUTPUT_PATH}/check_itsolver GMRES "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" 300 1e-7)
ADD_TEST(check_itsolver_rgmres ${EXECUTABLE_OUTPUT_PATH}/check_itsolver RGMRES "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" 50 1e-7)
ADD_TEST(check_itsolver_jacobi ${EXECUTABLE_OUTPUT_PATH}/check_itsolver JACOBI "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" 1000 1e-6)
ADD_TEST(check_itsolver_gaussseidel ${EXECUTABLE_OUTPUT_PATH}/check_itsolver GS "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" 1000 1e-6)
ADD_TEST(check_itsolver_sor ${EXECUTABLE_OUTPUT_PATH}/check_itsolver SOR "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" 1000 1e-6)
ADD_TEST(check_itsolver_ssor ${EXECUTABLE_OUTPUT_PATH}/check_itsolver SSOR "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" 1000 1e-6)

#Preconditioner
ADD_TEST(check_precond_diag ${EXECUTABLE_OUTPUT_PATH}/check_precond DIAG "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" 1000 1e-7)
ADD_TEST(check_precond_diag_cpx ${EXECUTABLE_OUTPUT_PATH}/check_precond DIAG "${CMAKE_SOURCE_DIR}/tests/data/randc50x50.mtx" 1000 1e-7)
ADD_TEST(check_precond_ilu0 ${EXECUTABLE_OUTPUT_PATH}/check_precond ILU0 "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" 1000 1e-7)
ADD_TEST(check_precond_iluk0 ${EXECUTABLE_OUTPUT_PATH}/check_precond ILUK0 "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" 1000 1e-7)
ADD_TEST(check_precond_iluk1 ${EXECUTABLE_OUTPUT_PATH}/check_precond ILUK1 "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" 1000 1e-7)
ADD_TEST(check_precond_iluk2 ${EXECUTABLE_OUTPUT_PATH}/check_precond ILUK2 "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" 1000 1e-7)
ADD_TEST(check_precond_iluk0_cpx ${EXECUTABLE_OUTPUT_PATH}/check_precond ILUK0 "${CMAKE_SOURCE_DIR}/tests/data/randc50x50.mtx" 1000 1e-7)
ADD_TEST(check_precond_iluk1_cpx ${EXECUTABLE_OUTPUT_PATH}/check_precond ILUK1 "${CMAKE_SOURCE_DIR}/tests/data/randc50x50.mtx" 1000 1e-7)
ADD_TEST(check_precond_iluk2_cpx ${EXECUTABLE_OUTPUT_PATH}/check_precond ILUK2 "${CMAKE_SOURCE_DIR}/tests/data/randc50x50.mtx" 1000 1e-7)
ADD_TEST(check_precond_sor ${EXECUTABLE_OUTPUT_PATH}/check_precond SOR "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" 1000 1e-7)
ADD_TEST(check_precond_ssor ${EXECUTABLE_OUTPUT_PATH}/check_precond SSOR "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" 1000 1e-7)






