###############################################################################
# set file locations
set(BASE_SRC_DIR base)
set(ALGOR_SRC_DIR algorithms)
set(FILECLASS_SRC_DIR fileclasses)
set(IMGCLASS_SRC_DIR imageclasses)
set(LASCLASS_SRC_DIR lasclasses)
set(APP_SRC_DIR apps)
set(PYTHON_SRC_DIR swig)

set(BASE_H
	${BASE_SRC_DIR}/IndexValue.h
	${BASE_SRC_DIR}/Optionpk.h
	${BASE_SRC_DIR}/PosValue.h
	${BASE_SRC_DIR}/Vector2d.h
	${BASE_SRC_DIR}/Vector2d.cc
	)

set(ALGOR_H
	${ALGOR_SRC_DIR}/ConfusionMatrix.h
	${ALGOR_SRC_DIR}/CostFactory.h
	${ALGOR_SRC_DIR}/CostFactorySVM.h
	${ALGOR_SRC_DIR}/Egcs.h
	${ALGOR_SRC_DIR}/FeatureSelector.h
	${ALGOR_SRC_DIR}/Filter.h
	${ALGOR_SRC_DIR}/Filter2d.h
	${ALGOR_SRC_DIR}/ImgRegression.h
	${ALGOR_SRC_DIR}/StatFactory.h
	${ALGOR_SRC_DIR}/myfann_cpp.h
	${ALGOR_SRC_DIR}/svm.h
	)

set(ALGOR_CC
	${ALGOR_SRC_DIR}/ConfusionMatrix.cc
	${ALGOR_SRC_DIR}/CostFactorySVM.cc
	${ALGOR_SRC_DIR}/CostFactorySVM.h
	${ALGOR_SRC_DIR}/Egcs.cc
	${ALGOR_SRC_DIR}/Filter.cc
	${ALGOR_SRC_DIR}/Filter2d.cc
	${ALGOR_SRC_DIR}/ImgRegression.cc
	${ALGOR_SRC_DIR}/svm.cpp
	)

set(FILECLASS_H
	${FILECLASS_SRC_DIR}/FileReaderAscii.h
	)

set(FILECLASS_CC
	${FILECLASS_SRC_DIR}/FileReaderAscii.cc
	)

set(IMGCLASS_H
	${IMGCLASS_SRC_DIR}/ImgRaster.h
	${IMGCLASS_SRC_DIR}/ImgCollection.h
	${IMGCLASS_SRC_DIR}/ImgReaderOgr.h
	${IMGCLASS_SRC_DIR}/ImgWriterOgr.h
	)

set(IMGCLASS_CC
	${IMGCLASS_SRC_DIR}/ImgRaster.cc
	${IMGCLASS_SRC_DIR}/ImgCollection.cc
	${IMGCLASS_SRC_DIR}/ImgReaderOgr.cc
	${IMGCLASS_SRC_DIR}/ImgWriterOgr.cc
	${IMGCLASS_SRC_DIR}/pkcomposite_lib.cc
	${IMGCLASS_SRC_DIR}/pkcrop_lib.cc
	${IMGCLASS_SRC_DIR}/pkextractogr_lib.cc
	${IMGCLASS_SRC_DIR}/pkextractimg_lib.cc
	${IMGCLASS_SRC_DIR}/pkstatprofile_lib.cc
	${IMGCLASS_SRC_DIR}/pkfilter_lib.cc
	${IMGCLASS_SRC_DIR}/pksvm_lib.cc
	${IMGCLASS_SRC_DIR}/pkdiff_lib.cc
	${IMGCLASS_SRC_DIR}/pkvalidate_lib.cc
	${IMGCLASS_SRC_DIR}/pkstretch_lib.cc
	)

set(LASCLASS_H
	${LASCLASS_SRC_DIR}/FileReaderLas.h
	)

set(LASCLASS_CC
	${LASCLASS_SRC_DIR}/FileReaderLas.cc
	)

set(APP_H
	${APP_SRC_DIR}/AppFactory.h
	)

set(PYTHON_CC
	${PYTHON_SRC_DIR}/ImgRaster_wrap.cc
	)

set(PKTOOLS_I
	${PYTHON_SRC_DIR}/pktools.i
	)
###############################################################################

###############################################################################
# Group source files for IDE source explorers
source_group("CMake Files" FILES CMakeLists.txt)
source_group("src_base" FILES ${BASE_H})
source_group("src_algor" FILES ${ALGOR_H} ${ALGOR_CC})
source_group("src_fileclass" FILES ${FILECLASS_H} ${FILECLASS_CC})
source_group("src_imgclass" FILES ${IMGCLASS_H} ${IMGCLASS_CC})
source_group("src_lasclass" FILES ${LASCLASS_H} ${LASCLASS_CC})
###############################################################################

###############################################################################
# Build and link library

add_library( ${PKTOOLS_BASE_LIB_NAME} ${BASE_H} )
target_link_libraries(${PKTOOLS_BASE_LIB_NAME} ${GDAL_LIBRARIES} ${GSL_LIBRARIES} )

add_library( ${PKTOOLS_IMAGECLASSES_LIB_NAME} ${IMGCLASS_H} ${IMGCLASS_CC} ${BASE_H} )
target_link_libraries(${PKTOOLS_IMAGECLASSES_LIB_NAME} ${GDAL_LIBRARIES} ${GSL_LIBRARIES} ${PKTOOLS_FILE_CLASSES_LIB_NAME})

add_library( ${PKTOOLS_ALGORITHMS_LIB_NAME} ${ALGOR_H} ${ALGOR_CC} ${FILECLASS_CC} ${FILECLASS_H} ${BASE_H} )
target_link_libraries(${PKTOOLS_ALGORITHMS_LIB_NAME} ${GDAL_LIBRARIES} ${GSL_LIBRARIES} ${PKTOOLS_IMAGECLASSES_LIB_NAME} ${PKTOOLS_FILE_CLASSES_LIB_NAME})

add_library( ${PKTOOLS_FILECLASSES_LIB_NAME} ${FILECLASS_H} ${FILECLASS_CC} ${BASE_H} )
target_link_libraries(${PKTOOLS_FILECLASSES_LIB_NAME} ${GDAL_LIBRARIES} ${GSL_LIBRARIES} )

if (BUILD_WITH_LIBLAS)
	add_library( ${PKTOOLS_LASCLASSES_LIB_NAME} ${LASCLASS_H} ${LASCLASS_CC} ${BASE_H} )
	target_link_libraries(${PKTOOLS_FILECLASSES_LIB_NAME} ${LIBLAS_LIBRARIES} ${BOOST_LIBRARIES} ${GDAL_LIBRARIES} ${GSL_LIBRARIES} )
endif(BUILD_WITH_LIBLAS)

###############################################################################
SET(PKTOOLS_LIBRARIES ${PKTOOLS_BASE_LIB_NAME} ${PKTOOLS_IMAGECLASSES_LIB_NAME} ${PKTOOLS_ALGORITHMS_LIB_NAME} ${PKTOOLS_FILECLASSES_LIB_NAME})
###############################################################################
# Set target properties
SET_TARGET_PROPERTIES(${PKTOOLS_BASE_LIB_NAME}
  PROPERTIES
  SOVERSION ${PKTOOLS_BASE_SOVERSION}
  VERSION ${PKTOOLS_BASE_VERSION}
)

SET_TARGET_PROPERTIES(${PKTOOLS_IMAGECLASSES_LIB_NAME}
  PROPERTIES
  SOVERSION ${PKTOOLS_IMAGECLASSES_SOVERSION}
  VERSION ${PKTOOLS_IMAGECLASSES_VERSION}
)

SET_TARGET_PROPERTIES(${PKTOOLS_ALGORITHMS_LIB_NAME}
  PROPERTIES
  SOVERSION ${PKTOOLS_ALGORITHMS_SOVERSION}
  VERSION ${PKTOOLS_ALGORITHMS_VERSION}
)

SET_TARGET_PROPERTIES(${PKTOOLS_FILECLASSES_LIB_NAME}
  PROPERTIES
  SOVERSION ${PKTOOLS_FILECLASSES_SOVERSION}
  VERSION ${PKTOOLS_FILECLASSES_VERSION}
)

if (BUILD_WITH_LIBLAS)
  SET_TARGET_PROPERTIES(${PKTOOLS_LASCLASSES_LIB_NAME}
    PROPERTIES
    SOVERSION ${PKTOOLS_LASCLASSES_SOVERSION}
    VERSION ${PKTOOLS_LASCLASSES_VERSION}
    )
  SET(PKTOOLS_LIBRARIES ${PKTOOLS_LIBRARIES} ${PKTOOLS_LASCLASSES_LIB_NAME})
endif(BUILD_WITH_LIBLAS)

  SET(PKTOOLS_LIBRARIES ${PKTOOLS_LIBRARIES} ${BASE_LIB})
###############################################################################

########## SWIG #############
IF(BUILD_WITH_PYTHON)
  FIND_PACKAGE(SWIG REQUIRED)
  FIND_PACKAGE(PythonLibs)
  # SET(PYTHON_SRC_DIR ${PROJECT_SOURCE_DIR}/swig)
  FILE(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/swig)
  INCLUDE(${SWIG_USE_FILE})
  INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR})
  INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
  SET(CMAKE_SWIG_FLAGS "")
  SET_PROPERTY(SOURCE ${PYTHON_SRC_DIR}/pktools.i PROPERTY CPLUSPLUS ON)
  SWIG_ADD_MODULE(pktools python ${PYTHON_SRC_DIR}/pktools.i ${IMGCLASS_H} ${IMGCLASS_CC})
  # SET(SWIG_MODULE_pktools_EXTRA_DEPS ${PKTOOLS_LIBRARIES})
  SET_SOURCE_FILES_PROPERTIES(${PYTHON_SRC_DIR}/pktools.i PROPERTIES SWIG_FLAGS "-Wall")
  SWIG_LINK_LIBRARIES(pktools ${PYTHON_LIBRARIES} ${PKTOOLS_LIBRARIES})
ENDIF(BUILD_WITH_PYTHON)

########## INSTALL ##########
INSTALL (TARGETS ${PKTOOLS_LIB_NAME} LIBRARY DESTINATION lib)
EXECUTE_PROCESS(COMMAND python -c "from distutils.sysconfig import get_python_lib; print get_python_lib()" OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
SET(PYTHON_INSTALL_DIR ${PYTHON_SITE_PACKAGES} CACHE PATH "Install directory for Python library")
INSTALL(TARGETS ${SWIG_MODULE_pktools_REAL_NAME} LIBRARY DESTINATION ${PYTHON_INSTALL_DIR})
INSTALL(FILES  ${CMAKE_CURRENT_BINARY_DIR}/pktools.py DESTINATION ${PYTHON_INSTALL_DIR})
INSTALL (FILES ${PYTHON_SRC_DIR}/pktools.i DESTINATION ${PROJECT_INCLUDE_DIR}/swig)
###############################################################################
###############################################################################
# Installation 
# install (TARGETS ${PKTOOLS_BASE_LIB_NAME} ${PKTOOLS_IMAGECLASSES_LIB_NAME} ${PKTOOLS_ALGORITHMS_LIB_NAME} ${PKTOOLS_FILECLASSES_LIB_NAME} DESTINATION lib)
INSTALL (TARGETS ${PKTOOLS_BASE_LIB_NAME} DESTINATION lib)
INSTALL (TARGETS ${PKTOOLS_IMAGECLASSES_LIB_NAME} DESTINATION lib)
INSTALL (TARGETS ${PKTOOLS_ALGORITHMS_LIB_NAME} DESTINATION lib)
INSTALL (TARGETS ${PKTOOLS_FILECLASSES_LIB_NAME} DESTINATION lib)
# INSTALL (TARGETS ${PKTOOLS_BASE_LIB_NAME} EXPORT base-targets DESTINATION lib)
# INSTALL (TARGETS ${PKTOOLS_IMAGECLASSES_LIB_NAME} EXPORT imageClasses-targets DESTINATION lib)
# INSTALL (TARGETS ${PKTOOLS_ALGORITHMS_LIB_NAME} EXPORT algorithms-targets DESTINATION lib)
# INSTALL (TARGETS ${PKTOOLS_FILECLASSES_LIB_NAME} EXPORT fileClasses-targets DESTINATION lib)

# INSTALL(EXPORT base-targets DESTINATION lib/cmake/pktools)
# INSTALL(EXPORT imageClasses-targets DESTINATION lib/cmake/pktools)
# INSTALL(EXPORT algorithms-targets DESTINATION lib/cmake/pktools)
# INSTALL(EXPORT fileClasses-targets DESTINATION lib/cmake/pktools)

INSTALL (FILES ${BASE_H} DESTINATION ${PROJECT_INCLUDE_DIR}/base)
INSTALL (FILES ${IMGCLASS_H} DESTINATION ${PROJECT_INCLUDE_DIR}/imageclasses)
INSTALL (FILES ${ALGOR_H} DESTINATION ${PROJECT_INCLUDE_DIR}/algorithms)
INSTALL (FILES ${FILECLASS_H} DESTINATION ${PROJECT_INCLUDE_DIR}/fileclasses)
INSTALL (FILES ${APP_H} DESTINATION ${PROJECT_INCLUDE_DIR}/apps)

if (BUILD_WITH_LIBLAS)
  INSTALL (TARGETS ${PKTOOLS_LASCLASSES_LIB_NAME} DESTINATION lib)
  INSTALL (FILES ${LASCLASS_H} DESTINATION ${PROJECT_INCLUDE_DIR}/lasclasses)
  # INSTALL(EXPORT lasclasses-targets DESTINATION lib/cmake/pktools)
endif(BUILD_WITH_LIBLAS)

###############################################################################
