cmake_minimum_required(VERSION 3.14 FATAL_ERROR)
project(BohdanAna)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)


find_package(ILCUTIL REQUIRED COMPONENTS ILCSOFT_CMAKE_MODULES)
include(ilcsoft_default_settings)

# Remove this ugly warning
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-effc++ -O3")

include_directories(${PROJECT_SOURCE_DIR}/include)
add_library(${PROJECT_NAME} SHARED ${PROJECT_SOURCE_DIR}/src/BohdanAna.cc
                                    ${PROJECT_SOURCE_DIR}/src/BohdanDrawing.cc
                                    ${PROJECT_SOURCE_DIR}/src/BohdanUtils.cc
                                    ${PROJECT_SOURCE_DIR}/src/TOF.cc
                                    ${PROJECT_SOURCE_DIR}/src/TrackLength.cc)

### DEPENDENCIES ###
find_package(Marlin REQUIRED)
find_package(MarlinUtil REQUIRED)
find_package(MarlinTrk REQUIRED)

FOREACH(pkg Marlin MarlinUtil MarlinTrk)
   include_directories(${${pkg}_INCLUDE_DIRS})
   target_link_libraries(${PROJECT_NAME} ${${pkg}_LIBRARIES})
ENDFOREACH()

# this is fixed in MarlinUtil 1.17.01
find_package( DD4hep REQUIRED COMPONENTS DDRec)
include_directories( SYSTEM ${DD4hep_INCLUDE_DIRS} )
target_link_libraries(${PROJECT_NAME} ${DD4hep_COMPONENT_LIBRARIES})

install(TARGETS ${PROJECT_NAME} DESTINATION ${PROJECT_SOURCE_DIR}/lib)
