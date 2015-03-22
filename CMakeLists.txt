cmake_minimum_required(VERSION 2.8)

project(Watershed)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

if (ITKVtkGlue_LOADED)
  find_package(VTK REQUIRED)
  include(${VTK_USE_FILE})
endif()

add_executable(WatershedExtractorTest watershed_extractor_test.cc
                                      watershed_extractor.cc
                                      util.cc)
target_link_libraries(WatershedExtractorTest ${VTK_LIBRARIES})

add_executable(ITKTest itk_test.cc)
target_link_libraries(ITKTest ${ITK_LIBRARIES})
target_link_libraries(ITKTest ${VTK_LIBRARIES})

add_executable(WatershedCorrectorTest watershed_corrector_test.cc
                                      watershed_corrector.cc
                                      util.cc)
target_link_libraries(WatershedCorrectorTest ${VTK_LIBRARIES})

