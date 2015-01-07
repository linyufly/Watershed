cmake_minimum_required(VERSION 2.8)

project(Watershed)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(WatershedExtractorTest watershed_extractor_test.cc
                                      watershed_extractor.cc
                                      util.cc)

target_link_libraries(WatershedExtractorTest ${VTK_LIBRARIES})
