cmake_minimum_required(VERSION 2.8.12.1)

PROJECT(graphic)
SET(VTK_DIR /Users/edward/VTK-8.1.2)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(graphic MACOSX_BUNDLE graphic)

if(VTK_LIBRARIES)
  target_link_libraries(graphic ${VTK_LIBRARIES})
else()
  target_link_libraries(graphic vtkHybrid)
endif()
