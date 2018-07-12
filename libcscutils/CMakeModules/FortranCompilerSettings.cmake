IF(NOT "${HOSTOPT}" STREQUAL "")
    IF(NOT (HOSTOPT STREQUAL OFF OR HOSTOPT STREQUAL ON))
        STRING(SUBSTRING "${HOSTOPT}" 0 1 FIRST_CHAR)
        IF( "${FIRST_CHAR}" STREQUAL "/" )
            MESSAGE(STATUS "Load user supplied Host Optimizations for Fortran -- ${HOSTOPT}")
            INCLUDE(${HOSTOPT})
        ELSE()
            MESSAGE(STATUS  "Load user supplied Host Optimizations for Fortran -- ${CMAKE_BINARY_DIR}/${HOSTOPT}")
            INCLUDE(${CMAKE_BINARY_DIR}/${HOSTOPT})
        ENDIF()
        SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${HOSTOPT_Fortran}")
    ENDIF()
ENDIF()



IF(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    # GNU
    SET(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -Wall -Wunused -Wimplicit-procedure")
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -frecursive -fPIC")

    IF(INTEGER8 STREQUAL ON)
        SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-integer-8")
    ENDIF()
    IF(HOSTOPT STREQUAL ON)
        SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -march=native -mtune=native")

    ENDIF()

    LIST(APPEND LIBRARIES "gfortran")

ELSEIF(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    # Intel
    SET(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -warn -g -warn nointerfaces")
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -recursive -fpic -fPIC -heap-arrays 256 ")

    IF(HOSTOPT STREQUAL ON)
        SET(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -O3 -xHost ")
        SET(CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS_DEBUG} -xHost")
    ENDIF()

    IF(INTEGER8 STREQUAL ON)
        SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -i8")
    ENDIF()
    LIST(APPEND LIBRARIES "ifcore")


ELSEIF(CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
    # PGI
    SET(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -O4 ")
    SET(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -g -Minfo=all")
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fPIC -fpic")
    IF(INTEGER8 STREQUAL ON)
        SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -i8")
    ENDIF()

ELSEIF(CMAKE_Fortran_COMPILER_ID STREQUAL "XL")
    # IBM XL
    IF(HOSTOPT STREQUAL ON)
        SET(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -O5 -qessl -qhot -qtune=auto -qarch=auto -qnoipa")
    ELSE()
        SET(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -O3 -qessl -qhot -qnoipa")
    ENDIF()
    SET(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -qstrict=ieeefp -qreport -qlistfmt=html=all -qnosave -qxlf77=nopersistent")
    SET(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -qstrict=ieeefp -qnosave -qxlf77=nopersistent")
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qpic -qnosave -qxlf77=nopersistent -qalias=std  -qnoipa")
    STRING(REPLACE "-qhalt=e" "" CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")

    IF(INTEGER8 STREQUAL ON)
        SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qintsize=8")
    ENDIF()

    IF(OPENMP_FOUND)
        LIST(REMOVE_ITEM CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES "xlomp_ser")
    ENDIF()
ENDIF()

