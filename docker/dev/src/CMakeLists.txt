# Copyright (c) General Electric Company, 2017.  All rights reserved.
#---------------------------------------------------------
# Wavelet Nucleii Segmentation Filter

#cmake_minimum_required(VERSION 3.4)

FIND_PACKAGE(ITK REQUIRED)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
  MESSAGE(STATUS "Using <${ITK_USE_FILE}>" )

	IF(ITK_SOURCE_DIR)
	  	MESSAGE(STATUS "ITK Source in <${ITK_SOURCE_DIR}>")
	ELSE(ITK_SOURCE_DIR)
		MESSAGE("Please set path to ITK_SOURCE_DIR")
  	 	FIND_PATH(ITK_SOURCE_DIR CMakeLists.txt)
 	ENDIF(ITK_SOURCE_DIR)  	
 	
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR
          "Cannot build without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)


ADD_EXECUTABLE( itkWaveletNucleiSegmentationTest 
		itkWaveletNucleiSegmentationTest.cxx
		itkWaveletNucleiSegmentationFilter.h
		itkWaveletNucleiSegmentationFilter.txx
		itkShapeWaterShedImageFilter.h
		itkShapeWaterShedImageFilter.txx
		 )

TARGET_LINK_LIBRARIES( 	itkWaveletNucleiSegmentationTest ${ITK_LIBRARIES})

INSTALL( TARGETS itkWaveletNucleiSegmentationTest
   RUNTIME DESTINATION bin
   LIBRARY DESTINATION lib
   ARCHIVE DESTINATION lib
)
