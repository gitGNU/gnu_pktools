###############################################################################
# set file locations
set(BASE_SRC_DIR base)
set(ALGOR_SRC_DIR algorithms)
set(FILECLASS_SRC_DIR fileclasses)
set(IMGCLASS_SRC_DIR imageclasses)
set(LASCLASS_SRC_DIR lasclasses)

set(BASE_H
	${BASE_SRC_DIR}/IndexValue.h
	${BASE_SRC_DIR}/Optionpk.h
	${BASE_SRC_DIR}/PosValue.h
	${BASE_SRC_DIR}/Vector2d.h
	${BASE_SRC_DIR}/Vector2d.cc
#	${BASE_SRC_DIR}/config.h
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
#	${ALGOR_SRC_DIR}/OptFactory.h
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
	${IMGCLASS_SRC_DIR}/ImgRasterGdal.h 
	${IMGCLASS_SRC_DIR}/ImgReaderGdal.h 
	${IMGCLASS_SRC_DIR}/ImgReaderOgr.h 
	${IMGCLASS_SRC_DIR}/ImgWriterGdal.h 
	${IMGCLASS_SRC_DIR}/ImgWriterOgr.h 
	${IMGCLASS_SRC_DIR}/ImgUpdaterGdal.h 
	)

set(IMGCLASS_CC
	${IMGCLASS_SRC_DIR}/ImgRasterGdal.cc
	${IMGCLASS_SRC_DIR}/ImgReaderGdal.cc
	${IMGCLASS_SRC_DIR}/ImgReaderOgr.cc
	${IMGCLASS_SRC_DIR}/ImgWriterGdal.cc 
	${IMGCLASS_SRC_DIR}/ImgWriterOgr.cc
	${IMGCLASS_SRC_DIR}/ImgUpdaterGdal.cc 
	)
	
set(LASCLASS_H
	${LASCLASS_SRC_DIR}/FileReaderLas.h 
	)

set(LASCLASS_CC
	${LASCLASS_SRC_DIR}/FileReaderLas.cc 
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
target_link_libraries(${PKTOOLS_BASE_LIB_NAME} ${GDAL_LIBRARIES} ${GSL_LIBRARIES} ${ARMADILLO_LIBRARIES} )

add_library( ${PKTOOLS_IMAGECLASSES_LIB_NAME} ${IMGCLASS_H} ${IMGCLASS_CC} ${BASE_H} )
target_link_libraries(${PKTOOLS_IMAGECLASSES_LIB_NAME} ${GDAL_LIBRARIES} ${GSL_LIBRARIES} ${ARMADILLO_LIBRARIES} )
# foreach(lib ${GSL_LIBRARIES})
# message(${lib})
# endforeach()
add_library( ${PKTOOLS_ALGORITHMS_LIB_NAME} ${ALGOR_H} ${ALGOR_CC} ${BASE_H} )
target_link_libraries(${PKTOOLS_ALGORITHMS_LIB_NAME} ${GDAL_LIBRARIES} ${GSL_LIBRARIES} ${ARMADILLO_LIBRARIES} ${PKTOOLS_IMAGECLASSES_LIB_NAME} )

add_library( ${PKTOOLS_FILECLASSES_LIB_NAME} ${FILECLASS_H} ${FILECLASS_CC} ${BASE_H} )
target_link_libraries(${PKTOOLS_FILECLASSES_LIB_NAME} ${GDAL_LIBRARIES} ${GSL_LIBRARIES} ${ARMADILLO_LIBRARIES} )

if (BUILD_WITH_LIBLAS)
	add_library( ${PKTOOLS_LASCLASSES_LIB_NAME} ${LASCLASS_H} ${LASCLASS_CC} ${BASE_H} )
	target_link_libraries(${PKTOOLS_FILECLASSES_LIB_NAME} ${LIBLAS_LIBRARIES} ${BOOST_LIBRARIES} ${GDAL_LIBRARIES} ${GSL_LIBRARIES} ${ARMADILLO_LIBRARIES} )
endif(BUILD_WITH_LIBLAS)

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
###############################################################################

###############################################################################
# Installation 
install (TARGETS ${PKTOOLS_BASE_LIB_NAME} ${PKTOOLS_IMAGECLASSES_LIB_NAME} ${PKTOOLS_ALGORITHMS_LIB_NAME} ${PKTOOLS_FILECLASSES_LIB_NAME} DESTINATION lib)

install (FILES ${BASE_H} DESTINATION ${PROJECT_INCLUDE_DIR}/base)
install (FILES ${IMGCLASS_H} DESTINATION ${PROJECT_INCLUDE_DIR}/imageclasses)
install (FILES ${ALGOR_H} DESTINATION ${PROJECT_INCLUDE_DIR}/algorithms)
install (FILES ${FILECLASS_H} DESTINATION ${PROJECT_INCLUDE_DIR}/fileclasses)

if (BUILD_WITH_LIBLAS)
install (TARGETS ${PKTOOLS_LASCLASSES_LIB_NAME} DESTINATION lib)
install (FILES ${LASCLASS_H} DESTINATION ${PROJECT_INCLUDE_DIR}/lasclasses)
endif(BUILD_WITH_LIBLAS)
###############################################################################
