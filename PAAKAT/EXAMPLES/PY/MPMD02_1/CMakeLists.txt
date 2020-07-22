cmake_minimum_required(VERSION 3.0)
project(PLEPP VERSION 0.1 LANGUAGES C CXX Fortran)  
set(CMAKE_CXX_STANDARD 11)

## VTK   
list(APPEND LIBS vtkPVCatalyst vtkPVVTKExtensionsDefault vtkPVClientServerCoreCore)
list(APPEND LIBS vtkFiltersParallelGeometry vtkFiltersParallelMPI vtkFiltersParallel)
list(APPEND LIBS vtkPVPythonCatalyst)
find_package(ParaView 4.1 REQUIRED COMPONENTS ${LIBS})
include("${PARAVIEW_USE_FILE}")


## PLEPP 
## Include Directories (for all targets)
set(PLEPP_ROOT "../../LIBPLEPP/" CACHE FILEPATH "path to src and include")
include_directories(${PLEPP_ROOT}/include) #${CMAKE_BINARY_DIR})

## Wildcard all the sources in src
file(GLOB PLEPP_SOURCES ${PLEPP_ROOT}/src/*.cpp ${PLEPP_ROOT}/src/*.c  *.cxx  *.f90)
add_library(PLEPP ${PLEPP_SOURCES})  

## Testing   
## 01 
set(exec01 iam_both.x)  
add_executable(${exec01} dummy.cxx)  
target_link_libraries(${exec01} PLEPP ${LIBS})   

## 03   
set(exec03 iam_alya.x)  
add_executable(${exec03} dummy.f90)  
set_property(TARGET ${exec03} PROPERTY LINKER_LANGUAGE Fortran)
target_link_libraries(${exec03} PLEPP ${LIBS})   

## 02   
set(exec02 iam_hemelb.x)
add_executable(${exec02} dummy.cxx)
target_compile_definitions(${exec02} PRIVATE -DWITHALYA)
target_link_libraries(${exec02} PLEPP ${LIBS}) 

