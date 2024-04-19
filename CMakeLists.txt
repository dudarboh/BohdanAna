cmake_minimum_required(VERSION 3.14 FATAL_ERROR)
project(BohdanAna LANGUAGES CXX)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

#Define target libraries
add_library(BohdanAna SHARED ${CMAKE_SOURCE_DIR}/src/BohdanAna.cc
                                    ${CMAKE_SOURCE_DIR}/src/BohdanDrawing.cc
                                    ${CMAKE_SOURCE_DIR}/src/BohdanUtils.cc
                                    ${CMAKE_SOURCE_DIR}/src/TOF.cc
                                    ${CMAKE_SOURCE_DIR}/src/TrackLength.cc)
add_library(CreateRefitPFO SHARED ${CMAKE_SOURCE_DIR}/src/CreateRefitPFO.cc)

### Dependencies ###
find_package(ILCUTIL REQUIRED COMPONENTS ILCSOFT_CMAKE_MODULES)
find_package( Marlin REQUIRED )
find_package( MarlinUtil REQUIRED )
find_package( MarlinTrk REQUIRED )
find_package( DD4hep REQUIRED COMPONENTS DDRec ) # Not needed if MarlinUtil >= 1.17.01
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-effc++ -O3")
include(ilcsoft_default_settings)

target_include_directories(BohdanAna PRIVATE ${Marlin_INCLUDE_DIRS}
${MarlinUtil_INCLUDE_DIRS}
${MarlinTrk_INCLUDE_DIRS}
${DD4hep_INCLUDE_DIRS}
${CMAKE_SOURCE_DIR}/include
)
target_link_libraries(BohdanAna PRIVATE ${Marlin_LIBRARIES}
                                       ${MarlinUtil_LIBRARIES}
                                       ${MarlinTrk_LIBRARIES}
                                       ${DD4hep_LIBRARIES}
)
target_compile_options(BohdanAna PRIVATE -Wno-effc++)

target_include_directories(CreateRefitPFO PRIVATE ${Marlin_INCLUDE_DIRS}
${MarlinUtil_INCLUDE_DIRS}
${CMAKE_SOURCE_DIR}/include
)
target_link_libraries(CreateRefitPFO PRIVATE ${Marlin_LIBRARIES}
                                       ${MarlinUtil_LIBRARIES}
)
target_compile_options(CreateRefitPFO PRIVATE -Wno-effc++)


install(TARGETS BohdanAna CreateRefitPFO DESTINATION ${CMAKE_SOURCE_DIR}/lib)
