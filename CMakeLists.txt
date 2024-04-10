CMAKE_MINIMUM_REQUIRED( VERSION 3.16 )
# policies not known to CMake 2.8 ...
cmake_policy(SET CMP0042 NEW) # Mac OS -rpath behavior

# declare the project name
PROJECT(genfit2)

# set project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 2 )
SET( ${PROJECT_NAME}_VERSION_MINOR 2 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )

include(GNUInstallDirs)

option(GF_BUILD_EVENT_DISPLAY "Build event display module." ON)

# install destinations can be passed via the command line:
# cmake -DLIBRARY_OUTPUT_DIRECTORY:PATH=.... -DINCLUDE_OUTPUT_DIRECTORY:PATH=....
# otherwise, set default install destinations
IF(NOT DEFINED LIBRARY_OUTPUT_DIRECTORY)
	set( LIBRARY_OUTPUT_DIRECTORY ${CMAKE_INSTALL_LIBDIR} )
ENDIF()
IF(NOT DEFINED INCLUDE_OUTPUT_DIRECTORY)
	set( INCLUDE_OUTPUT_DIRECTORY ${CMAKE_INSTALL_INCLUDEDIR} )
ENDIF()

# define path to macros and settings
LIST( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake )
# include the macros and settings
INCLUDE(settings)
INCLUDE(macros)
INCLUDE(genfit)

FIND_PACKAGE(ROOT REQUIRED COMPONENTS Geom)
IF(ROOT_VERSION VERSION_LESS 6.20.0)
	MESSAGE(FATAL_ERROR "GenFit requires a ROOT version > 6.20.00")
ENDIF()
set(CMAKE_CC_FLAGS "${CMAKE_CC_FLAGS} ${ROOT_CC_FLAGS}")


# IF(DEFINED ENV{RAVEPATH})
# 	MESSAGE(STATUS "Environment variable RAVEPATH is set. Try to build GFRave.")
# 	# get compiler flags from rave
# 	INCLUDE(FindPkgConfig)
# 	SET(ENV{PKG_CONFIG_PATH} $ENV{RAVEPATH})
# 	find_package(Rave REQUIRED)
	
# 	SET(RAVE True)
# ELSEIF(DEFINED Rave_LDFLAGS)
# 	IF(DEFINED Rave_INCLUDE_DIRS)
# 		IF(DEFINED Rave_CFLAGS)
# 			MESSAGE(STATUS "Rave_LDFLAGS, Rave_INCLUDE_DIRS and Rave_CFLAGS are set. Try to build GFRave.")
# 			SET(RAVE True)
# 		ENDIF()
# 	ENDIF()
# ELSE()
#    	MESSAGE(STATUS "No RAVEPATH or Rave flags given. GFRave will not be built.")
# ENDIF()

# IF(DEFINED RAVE)
# 	# turn comma separated lists into space separated strings
# 	string (REPLACE ";" " " Rave_LDFLAGS_STR "${Rave_LDFLAGS}")
# 	string (REPLACE ";" " " Rave_INCLUDE_DIRS_STR "${Rave_INCLUDE_DIRS}")
# 	string (REPLACE ";" " " Rave_CFLAGS_STR "${Rave_CFLAGS}")

# 	SET (CMAKE_CXX_FLAGS ${Rave_DEFINITIONS} ${Rave_CFLAGS_STR})
# 	SET(GF_INC_DIRS
# 		./GFRave/include/ ${Rave_INCLUDE_DIRS_STR}
# 	)
	
# 	AUX_SOURCE_DIRECTORY( ./GFRave/src  library_sources )
# ENDIF()

# Set GENFIT include directories
SET(GF_INC_DIRS
		./core/include/
		./fields/include/
		./finitePlanes/include/
		./fitters/include/
		./GBL/include/
		./measurements/include/
		./trackReps/include/
		./utilities/include/
		${GF_INC_DIRS}
)
if(GF_BUILD_EVENT_DISPLAY)
    list(APPEND GF_INC_DIRS ./eventDisplay/include/)
endif()

INSTALL(DIRECTORY ${GF_INC_DIRS}
        DESTINATION ${INCLUDE_OUTPUT_DIRECTORY}
        PATTERN ".svn" EXCLUDE)

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

# SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -Wpedantic -Wall -Wextra \
#  -Wno-long-long -Wshadow -Werror=overloaded-virtual -Wunused-parameter -Waddress")

# SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -O0 -g")
# SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3")
# SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS} -O3 -g")


INCLUDE_DIRECTORIES(
		BEFORE
		${ROOT_INCLUDE_DIRS}
		# SYSTEM ${EIGEN3_INCLUDE_DIR}
		${GF_INC_DIRS}
)

# FIND_PACKAGE(Doxygen)
# IF(DOXYGEN_FOUND)
#     CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile @ONLY)

#     ADD_CUSTOM_TARGET(doc
#         ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile
#         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc
#         COMMENT "Generating API documentation with Doxygen" VERBATIM)
# ELSE()
#     MESSAGE( STATUS "Doxygen not installed or found -- API documentation can't be generated!" )
# ENDIF(DOXYGEN_FOUND)

# add the list of source code files -- careful, when a new file is added!
AUX_SOURCE_DIRECTORY( ./core/src          library_sources )
if(GF_BUILD_EVENT_DISPLAY)
	AUX_SOURCE_DIRECTORY( ./eventDisplay/src  library_sources )
endif()
AUX_SOURCE_DIRECTORY( ./fields/src        library_sources )
AUX_SOURCE_DIRECTORY( ./finitePlanes/src  library_sources )
AUX_SOURCE_DIRECTORY( ./fitters/src       library_sources )
AUX_SOURCE_DIRECTORY( ./GBL/src           library_sources )
AUX_SOURCE_DIRECTORY( ./measurements/src  library_sources )
AUX_SOURCE_DIRECTORY( ./trackReps/src     library_sources )
AUX_SOURCE_DIRECTORY( ./utilities/src     library_sources )

# Dictionary generation.  For the time being, we list classes one-by-one.
SET(CORE_DICTIONARY_SOURCES
	AbsFinitePlane.h
	AbsFitterInfo.h
	AbsHMatrix.h
	AbsTrackRep.h
	DetPlane.h
	FitStatus.h
	IO.h
	Material.h
	MeasuredStateOnPlane.h
	MeasurementFactory.h
	MeasurementOnPlane.h
	MeasurementProducer.h
	StateOnPlane.h
	AbsMeasurement.h
	ThinScatterer.h
	Track.h
	TrackCand.h
	TrackCandHit.h
	TrackPoint.h
	RectangularFinitePlane.h)
ROOT_GENERATE_DICTIONARY(
	coreRootDict
	${CORE_DICTIONARY_SOURCES}
	LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/core/src/coreLinkDef.h)
SET(library_sources ${CMAKE_CURRENT_BINARY_DIR}/coreRootDict.cxx ${library_sources})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libcoreRootDict_rdict.pcm DESTINATION ${LIBRARY_OUTPUT_DIRECTORY})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libcoreRootDict.rootmap DESTINATION ${LIBRARY_OUTPUT_DIRECTORY})

if(GF_BUILD_EVENT_DISPLAY)
	SET(EVENTDISPLAY_DICTIONARY_SOURCES
		EventDisplay.h)
	ROOT_GENERATE_DICTIONARY(
		eventDisplayRootDict
		${EVENTDISPLAY_DICTIONARY_SOURCES}
		LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/eventDisplay/src/eventDisplayLinkDef.h)
	SET(library_sources ${CMAKE_CURRENT_BINARY_DIR}/eventDisplayRootDict.cxx ${library_sources})
	install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libeventDisplayRootDict_rdict.pcm DESTINATION ${LIBRARY_OUTPUT_DIRECTORY})
	install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libeventDisplayRootDict.rootmap DESTINATION ${LIBRARY_OUTPUT_DIRECTORY})
endif()

SET(FITTERS_DICTIONARY_SOURCES
	AbsKalmanFitter.h
	KalmanFittedStateOnPlane.h
	KalmanFitStatus.h
	ReferenceStateOnPlane.h
	KalmanFitterRefTrack.h
	KalmanFitter.h
	KalmanFitterInfo.h
	DAF.h)
ROOT_GENERATE_DICTIONARY(
	fittersRootDict
	${FITTERS_DICTIONARY_SOURCES}
	LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/fitters/src/fittersLinkDef.h)
SET(library_sources ${CMAKE_CURRENT_BINARY_DIR}/fittersRootDict.cxx ${library_sources})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libfittersRootDict_rdict.pcm DESTINATION ${LIBRARY_OUTPUT_DIRECTORY})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libfittersRootDict.rootmap DESTINATION ${LIBRARY_OUTPUT_DIRECTORY})

SET(GBL_DICTIONARY_SOURCES
	GFGbl.h
	GblFitter.h
	ICalibrationParametersDerivatives.h
	GblFitterInfo.h
	GblFitStatus.h
	GblData.h)
ROOT_GENERATE_DICTIONARY(
	GBLRootDict
    ${GBL_DICTIONARY_SOURCES}
	LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/GBL/src/GBLLinkDef.h)
SET(library_sources ${CMAKE_CURRENT_BINARY_DIR}/GBLRootDict.cxx ${library_sources})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libGBLRootDict_rdict.pcm DESTINATION ${LIBRARY_OUTPUT_DIRECTORY})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libGBLRootDict.rootmap DESTINATION ${LIBRARY_OUTPUT_DIRECTORY})

SET(MEASUREMENTS_DICTIONARY_SOURCES
	HMatrixPhi.h
	HMatrixU.h
	HMatrixUnit.h
	HMatrixV.h
	HMatrixUV.h
	FullMeasurement.h
	PlanarMeasurement.h
	ProlateSpacepointMeasurement.h
	SpacepointMeasurement.h
	WireMeasurement.h
	WireMeasurementNew.h
	WirePointMeasurement.h
	WireTrackCandHit.h)
ROOT_GENERATE_DICTIONARY(
	measurementsRootDict
	${MEASUREMENTS_DICTIONARY_SOURCES}
	LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/measurements/src/measurementsLinkDef.h)
SET(library_sources ${CMAKE_CURRENT_BINARY_DIR}/measurementsRootDict.cxx ${library_sources})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libmeasurementsRootDict_rdict.pcm DESTINATION ${LIBRARY_OUTPUT_DIRECTORY})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libmeasurementsRootDict.rootmap DESTINATION ${LIBRARY_OUTPUT_DIRECTORY})

SET(TRACKREPS_DICTIONARY_SOURCES
	RKTrackRep.h
	MplTrackRep.h)
ROOT_GENERATE_DICTIONARY(
	trackRepsRootDict
	${TRACKREPS_DICTIONARY_SOURCES}
	LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/trackReps/src/trackRepsLinkDef.h)
SET(library_sources ${CMAKE_CURRENT_BINARY_DIR}/trackRepsRootDict.cxx ${library_sources})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libtrackRepsRootDict_rdict.pcm DESTINATION ${LIBRARY_OUTPUT_DIRECTORY})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libtrackRepsRootDict.rootmap DESTINATION ${LIBRARY_OUTPUT_DIRECTORY})

SET(UTILITIES_DICTIONARY_SOURCES
	HelixTrackModel.h
	MeasurementCreator.h
	mySpacepointDetectorHit.h
	mySpacepointMeasurement.h)
ROOT_GENERATE_DICTIONARY(
	utilitiesRootDict
	${UTILITIES_DICTIONARY_SOURCES}
	LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/utilities/src/utilitiesLinkDef.h)
SET(library_sources ${CMAKE_CURRENT_BINARY_DIR}/utilitiesRootDict.cxx ${library_sources})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libutilitiesRootDict_rdict.pcm DESTINATION ${LIBRARY_OUTPUT_DIRECTORY})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libutilitiesRootDict.rootmap DESTINATION ${LIBRARY_OUTPUT_DIRECTORY})

if(DEFINED RAVE)
	SET(GFRAVE_DICTIONARY_SOURCES
		GFRaveTrackParameters.h
		GFRaveVertex.h)
	ROOT_GENERATE_DICTIONARY(
		GFRaveRootDict
		${GFRAVE_DICTIONARY_SOURCES}
		LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/GFRave/src/GFRaveLinkDef.h)
	SET(library_sources ${CMAKE_CURRENT_BINARY_DIR}/GFRaveRootDict.cxx ${library_sources})
	install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libGFRaveRootDict_rdict.pcm DESTINATION ${LIBRARY_OUTPUT_DIRECTORY})
	install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libGFRaveRootDict.rootmap DESTINATION ${LIBRARY_OUTPUT_DIRECTORY})
endif()


# install
ADD_SHARED_LIBRARY( ${PROJECT_NAME} ${library_sources})
INSTALL(TARGETS ${PROJECT_NAME} LIBRARY DESTINATION ${LIBRARY_OUTPUT_DIRECTORY})

# link ROOT libraries
if(GF_BUILD_EVENT_DISPLAY)
	set(GF_LINK_ROOT_LIBRARIES ROOT::Core ROOT::Physics ROOT::EG ROOT::Eve)
else()
	set(GF_LINK_ROOT_LIBRARIES ROOT::Core ROOT::Physics ROOT::EG)
endif()
if(DEFINED RAVE)
	TARGET_LINK_LIBRARIES(
			${PROJECT_NAME}
			${GF_LINK_ROOT_LIBRARIES}
			${Rave_LIB}
			${Rave_LDFLAGS_STR}
	)
else()
	TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${GF_LINK_ROOT_LIBRARIES})
endif()
unset(GF_LINK_ROOT_LIBRARIES)

# INCLUDE(CTest)

# ADD_CUSTOM_TARGET( tests )

# ADD_GENFIT_TEST( fitterTests               test/fitterTests/main.cc)
# ADD_GENFIT_TEST( minimalFittingExample     test/minimalFittingExample/main.cc)
# ADD_GENFIT_TEST( minimalFittingExample2    test/minimalFittingExample/main2.cc)
# ADD_GENFIT_TEST( measurementFactoryExample test/measurementFactoryExample/main.cc)
# ADD_GENFIT_TEST( streamerTest              test/streamerTest/main.cc)
# ADD_GENFIT_TEST( unitTests                 test/unitTests/main.cc)
# IF(DEFINED RAVE)
#   ADD_GENFIT_TEST( vertexingTest           test/vertexingTest/main.cc)
#   ADD_GENFIT_TEST( vertexingTestRead       test/vertexingTest/read.cc)
# ENDIF()

# IF(BUILD_TESTING)

# 	# Run these automatically via ctest since they do not need a GUI. Create the
# 	# necessary geometry first via a dummy test
# 	ADD_TEST(NAME createGeometry COMMAND root -l -b -q ${CMAKE_CURRENT_LIST_DIR}/test/makeGeom.C)

# 	ADD_TEST(NAME streamerTest COMMAND streamerTest)
# 	SET_PROPERTY(TEST streamerTest PROPERTY DEPENDS createGeometry)

# 	ADD_TEST(NAME unitTests COMMAND unitTests)
# 	SET_PROPERTY(TEST unitTests PROPERTY DEPENDS createGeometry)

#     ENABLE_TESTING()


# 	FIND_PACKAGE(GTest REQUIRED)
# 	ADD_EXECUTABLE(gtests
# 			gtest/TestConstField.cpp
# 			gtest/TestMaterialEffects.cpp
# 			gtest/TestRKTrackRep.cpp
# 			gtest/TestMaterial.cpp
# 			)
# 	TARGET_LINK_LIBRARIES(gtests ${GTEST_BOTH_LIBRARIES} ${ROOT_LIBS} ${PROJECT_NAME})  # gtest gtest_main
# 	MESSAGE(STATUS  ${GTEST_INCLUDE_DIRS})
# 	TARGET_INCLUDE_DIRECTORIES(gtests PUBLIC ${GTEST_INCLUDE_DIRS})
# 	GTEST_DISCOVER_TESTS(gtests)
# 	INSTALL(TARGETS gtests DESTINATION bin)
# ENDIF()

# generate and install following configuration files
#GENERATE_PACKAGE_CONFIGURATION_FILES( genfitConfig.cmake genfitConfigVersion.cmake )
