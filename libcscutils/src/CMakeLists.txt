SET(GENERIC_LIB_VERSION "1.0")
SET(GENERIC_LIB_SOVERSION "1")


SET(COMMON_SRC
    common/error_message.c
    common/backtrace.c
    common/counter.c
    common/string.c
    common/mtx_write.c
    common/table.c
    common/table_comment.c
    common/table_latex.c
    common/fortran_api.f90
    common/table.f90
    common/statistics.f90
)

IF(IO_ENABLE)
    SET(IO_SRC
        file/compress_handler.c
        file/compress_io_bzip.c
        file/compress_io_gzip.c
        file/compress_io_uncompressed.c
        file/compress_io_xz.c
        file/io.c
        file/dir.c
        )
ENDIF()

IF(INIFILE_ENABLE)
    SET(INI_SRC
        ini/inifile.c)
ENDIF()

IF(THREADING_ENABLE)
    SET(THREAD_SRC
        threading/worker_queue.c
        threading/bgthread.c)
ENDIF()

IF(HDF5_ENABLE)
    SET(HDF5_SRC
        hdf5/common.c
        hdf5/vector.c
        hdf5/matrix.c
        hdf5/datatypes.c
        hdf5/sparse_matrix.c
        hdf5/fortran_api.f90
        )
ENDIF()

IF(IMAGE_ENABLE)
    SET(IMAGE_SRC
        image/bmp.c)
ENDIF()

IF(HARDWARE_ENABLE)
    SET(HARDWARE_SRC
        hardware/meminfo.c
        hardware/cpucount.c
        hardware/cpufreq.c)
ENDIF()

IF( DS_ENABLE )
    SET(DS_SRC
        ds/ds.c
        ds/hashtable.c
        ds/slist.c
        )
ENDIF()

IF(TRACE_ENABLE)
    SET(TRACE_SRC
        eventtracer/reset.c
        eventtracer/write.c
        eventtracer/record_event.c
        eventtracer/fortran_api.f90
    )
ENDIF()

SET(SRC ${COMMON_SRC} ${FILE_SRC} ${IO_SRC} ${INI_SRC} ${THREAD_SRC} ${HDF5_SRC} ${IMAGE_SRC} ${HARDWARE_SRC} ${DS_SRC} ${TRACE_SRC})

INCLUDE_DIRECTORIES(${INCLUDE_DIR})
#INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR}/include)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/src)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/include)


ADD_LIBRARY(cscutils STATIC ${SRC})
TARGET_LINK_LIBRARIES(cscutils  ${LIBS})
SET_TARGET_PROPERTIES(cscutils PROPERTIES OUTPUT_NAME "cscutils")

# Build Shared Library if we do not use a stand alone build
IF(NOT INPROJECT)
    ADD_LIBRARY(cscutils_shared SHARED ${SRC})
    TARGET_LINK_LIBRARIES(cscutils_shared  ${LIBS})
    SET_TARGET_PROPERTIES(cscutils_shared  PROPERTIES VERSION ${GENERIC_LIB_VERSION} SOVERSION ${GENERIC_LIB_SOVERSION} OUTPUT_NAME "cscutils")
ENDIF()




IF(NOT INPROJECT)
    INSTALL(TARGETS cscutils cscutils_shared
        RUNTIME DESTINATION bin COMPONENT bin
        LIBRARY DESTINATION lib COMPONENT shlib
        ARCHIVE DESTINATION lib COMPONENT stlib
        PUBLIC_HEADER DESTINATION include COMPONENT dev)

    INSTALL(DIRECTORY ${PROJECT_SOURCE_DIR}/include/cscutils/ DESTINATION include/cscutils
        FILES_MATCHING PATTERN "*.h")
ENDIF()
