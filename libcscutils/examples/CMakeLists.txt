SET(SRC
    demo_inifile
    demo_io
    demo_table
    )
SET(FSRC
    demo_table_fortran
)

IF(HDF5_ENABLE)
    INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIRS})
    INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR}/src)
    SET(SRC ${SRC} demo_hdf)
    SET(FSRC ${FSRC} demo_hdf_f)


ENDIF()
IF(HARDWARE_ENABLE)
    SET(SRC ${SRC} demo_hardware)
ENDIF()

IF(DS_ENABLE)
    SET(SRC ${SRC} demo_hashtable
        demo_slist)

ENDIF()
IF(THREADING_ENABLE)
    SET(SRC ${SRC}
        demo_background_worker)
ENDIF()


INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/include)
FOREACH(srcfile ${SRC})
    ADD_EXECUTABLE(${srcfile} ${srcfile}.c)
    TARGET_LINK_LIBRARIES(${srcfile} cscutils ${LIBS})
ENDFOREACH()

FOREACH(srcfile ${FSRC})
    ADD_EXECUTABLE(${srcfile} ${srcfile}.f90)
    TARGET_LINK_LIBRARIES(${srcfile} cscutils ${LIBS})
ENDFOREACH()


