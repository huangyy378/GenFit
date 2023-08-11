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
IF(${ROOT_VERSION} VERSION_LESS 6.00.00)
	MESSAGE(FATAL_ERROR "GenFit requires a ROOT version > 6.00.00")
ENDIF()

# Set GENFIT include directories
SET(GF_INC_DIRS
		./core/include/
		./eventDisplay/include/
		./fields/include/
		./finitePlanes/include/
		./fitters/include/
		./GBL/include/
		./measurements/include/
		./trackReps/include/
		./utilities/include/
		${GF_INC_DIRS}
)

INSTALL(DIRECTORY ${GF_INC_DIRS}
        DESTINATION ${INCLUDE_OUTPUT_DIRECTORY}
        PATTERN ".svn" EXCLUDE)

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

INCLUDE_DIRECTORIES(
		BEFORE
		${ROOT_INCLUDE_DIRS}
		${GF_INC_DIRS}
)

# declare the ROOT library path
LINK_DIRECTORIES(${ROOT_LIBRARY_DIR})

# add the list of source code files -- careful, when a new file is added!
AUX_SOURCE_DIRECTORY( ./core/src          library_sources )
AUX_SOURCE_DIRECTORY( ./eventDisplay/src  library_sources )
AUX_SOURCE_DIRECTORY( ./fields/src        library_sources )
AUX_SOURCE_DIRECTORY( ./finitePlanes/src  library_sources )
AUX_SOURCE_DIRECTORY( ./fitters/src       library_sources )
AUX_SOURCE_DIRECTORY( ./GBL/src           library_sources )
AUX_SOURCE_DIRECTORY( ./measurements/src  library_sources )
AUX_SOURCE_DIRECTORY( ./trackReps/src     library_sources )
AUX_SOURCE_DIRECTORY( ./utilities/src     library_sources )

# add target
ADD_SHARED_LIBRARY( ${PROJECT_NAME} ${library_sources})

# Dictionary generation.  For the time being, we list classes one-by-one.
SET(DICTIONARY_SOURCES
		${CMAKE_CURRENT_SOURCE_DIR}/core/include/AbsFinitePlane.h
		${CMAKE_CURRENT_SOURCE_DIR}/core/include/AbsFitterInfo.h
		${CMAKE_CURRENT_SOURCE_DIR}/core/include/AbsHMatrix.h
		${CMAKE_CURRENT_SOURCE_DIR}/core/include/AbsTrackRep.h
		${CMAKE_CURRENT_SOURCE_DIR}/core/include/DetPlane.h
		${CMAKE_CURRENT_SOURCE_DIR}/core/include/FitStatus.h
		${CMAKE_CURRENT_SOURCE_DIR}/core/include/IO.h
		${CMAKE_CURRENT_SOURCE_DIR}/core/include/Material.h
		${CMAKE_CURRENT_SOURCE_DIR}/core/include/MeasuredStateOnPlane.h
		${CMAKE_CURRENT_SOURCE_DIR}/core/include/MeasurementFactory.h
		${CMAKE_CURRENT_SOURCE_DIR}/core/include/MeasurementOnPlane.h
		${CMAKE_CURRENT_SOURCE_DIR}/core/include/MeasurementProducer.h
		${CMAKE_CURRENT_SOURCE_DIR}/core/include/StateOnPlane.h
		${CMAKE_CURRENT_SOURCE_DIR}/core/include/AbsMeasurement.h
		${CMAKE_CURRENT_SOURCE_DIR}/core/include/ThinScatterer.h
		${CMAKE_CURRENT_SOURCE_DIR}/core/include/Track.h
		${CMAKE_CURRENT_SOURCE_DIR}/core/include/TrackCand.h
		${CMAKE_CURRENT_SOURCE_DIR}/core/include/TrackCandHit.h
		${CMAKE_CURRENT_SOURCE_DIR}/core/include/TrackPoint.h
		${CMAKE_CURRENT_SOURCE_DIR}/finitePlanes/include/RectangularFinitePlane.h
		
		${CMAKE_CURRENT_SOURCE_DIR}/eventDisplay/include/EventDisplay.h

		${CMAKE_CURRENT_SOURCE_DIR}/fitters/include/AbsKalmanFitter.h
		${CMAKE_CURRENT_SOURCE_DIR}/fitters/include/KalmanFittedStateOnPlane.h
		${CMAKE_CURRENT_SOURCE_DIR}/fitters/include/KalmanFitStatus.h
		${CMAKE_CURRENT_SOURCE_DIR}/fitters/include/ReferenceStateOnPlane.h
		${CMAKE_CURRENT_SOURCE_DIR}/fitters/include/KalmanFitterRefTrack.h
		${CMAKE_CURRENT_SOURCE_DIR}/fitters/include/KalmanFitter.h
		${CMAKE_CURRENT_SOURCE_DIR}/fitters/include/KalmanFitterInfo.h
		${CMAKE_CURRENT_SOURCE_DIR}/fitters/include/DAF.h

		${CMAKE_CURRENT_SOURCE_DIR}/GBL/include/GFGbl.h
		${CMAKE_CURRENT_SOURCE_DIR}/GBL/include/GblFitter.h
		${CMAKE_CURRENT_SOURCE_DIR}/GBL/include/ICalibrationParametersDerivatives.h
		${CMAKE_CURRENT_SOURCE_DIR}/GBL/include/GblFitterInfo.h
		${CMAKE_CURRENT_SOURCE_DIR}/GBL/include/GblFitStatus.h
		${CMAKE_CURRENT_SOURCE_DIR}/GBL/include/GblData.h

		${CMAKE_CURRENT_SOURCE_DIR}/measurements/include/HMatrixPhi.h
		${CMAKE_CURRENT_SOURCE_DIR}/measurements/include/HMatrixU.h
		${CMAKE_CURRENT_SOURCE_DIR}/measurements/include/HMatrixUnit.h
		${CMAKE_CURRENT_SOURCE_DIR}/measurements/include/HMatrixV.h
		${CMAKE_CURRENT_SOURCE_DIR}/measurements/include/HMatrixUV.h
		${CMAKE_CURRENT_SOURCE_DIR}/measurements/include/FullMeasurement.h
		${CMAKE_CURRENT_SOURCE_DIR}/measurements/include/PlanarMeasurement.h
		${CMAKE_CURRENT_SOURCE_DIR}/measurements/include/ProlateSpacepointMeasurement.h
		${CMAKE_CURRENT_SOURCE_DIR}/measurements/include/SpacepointMeasurement.h
		${CMAKE_CURRENT_SOURCE_DIR}/measurements/include/WireMeasurement.h
		${CMAKE_CURRENT_SOURCE_DIR}/measurements/include/WireMeasurementNew.h
		${CMAKE_CURRENT_SOURCE_DIR}/measurements/include/WirePointMeasurement.h
		${CMAKE_CURRENT_SOURCE_DIR}/measurements/include/WireTrackCandHit.h

		${CMAKE_CURRENT_SOURCE_DIR}/trackReps/include/RKTrackRep.h
		${CMAKE_CURRENT_SOURCE_DIR}/trackReps/include/MplTrackRep.h

		${CMAKE_CURRENT_SOURCE_DIR}/utilities/include/HelixTrackModel.h
		${CMAKE_CURRENT_SOURCE_DIR}/utilities/include/MeasurementCreator.h
		${CMAKE_CURRENT_SOURCE_DIR}/utilities/include/mySpacepointDetectorHit.h
		${CMAKE_CURRENT_SOURCE_DIR}/utilities/include/mySpacepointMeasurement.h
)
ROOT_GENERATE_DICTIONARY(
		${PROJECT_NAME}RootDict
		${CORE_DICTIONARY_SOURCES}
		LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/core/src/coreLinkDef.h 
		MODULE ${PROJECT_NAME}
)

# install
INSTALL(TARGETS ${PROJECT_NAME} LIBRARY DESTINATION ${LIBRARY_OUTPUT_DIRECTORY})

# link ROOT libraries
TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${ROOT_LIBRARIES})
