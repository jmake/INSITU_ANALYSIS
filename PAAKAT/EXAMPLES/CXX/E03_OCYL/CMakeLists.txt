cmake_minimum_required(VERSION 3.3)
project(nek5k01_01 CXX C Fortran)

## VTK   
list(APPEND LIBS vtkPVCatalyst vtkPVVTKExtensionsDefault vtkPVClientServerCoreCore)
list(APPEND LIBS vtkFiltersParallelGeometry vtkFiltersParallelMPI vtkFiltersParallel)
#list(APPEND LIBS vtkPVPythonCatalyst) # Only for python 
#set(CMAKE_CXX_FLAGS "-DWITH_PYTHON")  # Only for python 
find_package(ParaView 4.1 REQUIRED COMPONENTS ${LIBS})
include("${PARAVIEW_USE_FILE}")

## Testing   
FILE(GLOB SRCs *.cxx)
add_library(PAAKAT ${SRCs})


## Testing   
set(exec01 nek5k01_01.x)
add_executable(${exec01} nek5k01_01.f90)  
set_property(TARGET ${exec01} PROPERTY LINKER_LANGUAGE Fortran)
target_link_libraries(${exec01} PAAKAT ${LIBS})

