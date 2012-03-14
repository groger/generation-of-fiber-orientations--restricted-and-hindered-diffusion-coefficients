cmake_minimum_required(VERSION 2.6)

project(FiberSig)

option(BUILD_SHARED_LIBS "Build with shared libraries." OFF)

option(BUILD_FIBERODF_AS_LIBRARY "Build fiberodf as a library rather than an executable." OFF)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

find_package(Teem REQUIRED)
include(${Teem_USE_FILE})

find_package(VXL REQUIRED)
include(${VXL_CMAKE_DIR}/UseVXL.cmake)

add_definitions(-Wall)

include_directories(${PROJECT_SOURCE_DIR}/Common)
include_directories(${PROJECT_SOURCE_DIR}/Geometry)
include_directories(${PROJECT_SOURCE_DIR})
add_subdirectory(Common)
add_subdirectory(Geometry)

FIND_PACKAGE(GenerateCLP REQUIRED)
IF(GenerateCLP_FOUND)
  INCLUDE(${GenerateCLP_USE_FILE})
ELSE(GenerateCLP_FOUND)
  MESSAGE(FATAL_ERROR
          "GenerateCLP not found. Please set GenerateCLP_DIR.")
ENDIF(GenerateCLP_FOUND)

set(FIBERODF_SOURCE SphereIkosahedron.h Counter.h Counter.cc CounterStatistics.h CounterStatistics.cc Cylinder.h Cylinder.cc FiberODF_Common.h FiberODF_Common.cc FiberSegment.h FiberSegment.cc FiberSegmentExtractor.h FiberSegmentExtractor.cc GaussianOnSphere.h GaussianOnSphere.cc VolumeCounter.h VolumeCounter.cc DisplayFunc.h Stats.h FiberOrientation.h)
GENERATECLP(FIBERODF_SOURCE main.xml)

if(BUILD_FIBERODF_AS_LIBRARY)
  add_library(fiberodf ${FIBERODF_SOURCE})
else(BUILD_FIBERODF_AS_LIBRARY)
  add_executable(fibersig ${FIBERODF_SOURCE} main.cc)
endif(BUILD_FIBERODF_AS_LIBRARY)

target_link_libraries(fibersig common geometry vtkIO vtkGraphics ${ITK_LIBRARIES} teem vnl)
