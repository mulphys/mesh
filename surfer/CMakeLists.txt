
cmake_minimum_required(VERSION 2.8)

PROJECT(Surfer)
enable_language (Fortran)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

#set (CMAKE_Fortran_FLAGS_RELEASE "-traditional-cpp -funroll-all-loops -fno-f2c -O3")
#set (CMAKE_Fortran_FLAGS_DEBUG   "-cpp -fno-f2c -O0 -g")
#set (CMAKE_Fortran_FLAGS_DEBUG   "-traditional-cpp -O0 -g -fbacktrace -ffpe-trap=zero,overflow,underflow")
#set (CMAKE_Fortran_FLAGS_DEBUG   "-traditional-cpp -O0 -g -fbacktrace -ffpe-trap=list")
#set (CMAKE_Fortran_FLAGS_DEBUG   "-traditional-cpp")
#SET(CMAKE_CXX_FLAGS_DEBUG  "-O0 -g")
#SET(CMAKE_C_FLAGS_DEBUG  "-O0 -g")

add_library(lib_io SHARED io.c)
add_library(lib_vtk SHARED VTK.cc)
add_library(lib_fsurfer SHARED surfer.f)

if(VTK_LIBRARIES)
  target_link_libraries(lib_vtk ${VTK_LIBRARIES})
else()
  target_link_libraries(lib_vtk vtkHybrid vtkWidgets)
endif()

add_executable(surfer surfer.f90)
target_link_libraries(surfer lib_fsurfer lib_io lib_vtk)

add_executable(frontmesh frontmesh.cc frontmesh.h)
# add this if using VTK output
#target_link_libraries(frontmesh lib_vtk)

add_executable(tetgen2stl MACOSX_BUNDLE tetgen2stl.cc)
target_link_libraries(tetgen2stl lib_vtk)
