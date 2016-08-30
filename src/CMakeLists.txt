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

# if (BUILD_WITH_PYTHON)
#   FIND_PACKAGE(SWIG REQUIRED)
#   INCLUDE(${SWIG_USE_FILE})  
#   FIND_PACKAGE(PythonLibs)
#   INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
#   INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR} )
#   INCLUDE_DIRECTORIES(${IMGCLASS_SRC_DIR})
#   SET(CMAKE_SWIG_FLAGS "")

#   SET_SOURCE_FILES_PROPERTIES(ImgRaster.i PROPERTIES CPLUSPLUS ON)
#   SWIG_ADD_MODULE(ImgRaster python ${PYTHON_I} ${PYTHON_CC})
#   SWIG_LINK_LIBRARIES(ImgRaster ${PKTOOLS_IMAGECLASSES_LIB_NAME} ${PYTHON_LIBRARIES})
# #  SET_SOURCE_FILES_PROPERTIES(${PYTHON_I} PROPERTIES SWIG_FLAGS "-includeall")

# #  message("add module")
# #  SWIG_ADD_MODULE(${PKTOOLS_PYTHON_MODULE} python ${PYTHON_I})
# #  SWIG_LINK_LIBRARIES(${PKTOOLS_PYTHON_MODULE} ${PKTOOLS_IMAGECLASSES_LIB_NAME} ${PYTHON_LIBRARIES})
# endif(BUILD_WITH_PYTHON)
###############################################################################

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
endif(BUILD_WITH_LIBLAS)

# if (BUILD_WITH_PYTHON)
#   SET_TARGET_PROPERTIES(${PKTOOLS_PYTHON_MODULE}
#     PROPERTIES
#     SOVERSION ${PKTOOLS_PYTHON_SOVERSION}
#     VERSION ${PKTOOLS_PYTHON_VERSION}
#     )
# endif(BUILD_WITH_PYTHON)

###############################################################################

###############################################################################
# Installation 
install (TARGETS ${PKTOOLS_BASE_LIB_NAME} ${PKTOOLS_IMAGECLASSES_LIB_NAME} ${PKTOOLS_ALGORITHMS_LIB_NAME} ${PKTOOLS_FILECLASSES_LIB_NAME} DESTINATION lib)

install (FILES ${BASE_H} DESTINATION ${PROJECT_INCLUDE_DIR}/base)
install (FILES ${IMGCLASS_H} DESTINATION ${PROJECT_INCLUDE_DIR}/imageclasses)
install (FILES ${ALGOR_H} DESTINATION ${PROJECT_INCLUDE_DIR}/algorithms)
install (FILES ${FILECLASS_H} DESTINATION ${PROJECT_INCLUDE_DIR}/fileclasses)
install (FILES ${APP_H} DESTINATION ${PROJECT_INCLUDE_DIR}/apps)

if (BUILD_WITH_LIBLAS)
  install (TARGETS ${PKTOOLS_LASCLASSES_LIB_NAME} DESTINATION lib)
  install (FILES ${LASCLASS_H} DESTINATION ${PROJECT_INCLUDE_DIR}/lasclasses)
endif(BUILD_WITH_LIBLAS)

if (BUILD_WITH_PYTHON)
   install (FILES ${PKTOOLS_I} DESTINATION ${PROJECT_INCLUDE_DIR}/swig)
#   #  SET_SOURCE_FILES_PROPERTIES(ImgRaster.i PROPERTIES SWIG_FLAGS "-includeall")
#   SWIG_ADD_MODULE(ImgRaster python ${PYTHON_CC} ${PYTHON_I})
#   SWIG_LINK_LIBRARIES(ImgRaster ${PKTOOLS_IMAGECLASSES_LIB_NAME} ${PYTHON_LIBRARIES})
#   #INSTALL(TARGETS ${PKTOOLS_PYTHON_MODULE} DESTINATION ${PYTHON_SITE_PACKAGES})
endif(BUILD_WITH_PYTHON)
###############################################################################
