cmake_minimum_required(VERSION 3.0) 
project(PLEPP VERSION 0.1 LANGUAGES C CXX)  
set(CMAKE_CXX_STANDARD 11)

## 1.1. SOURCES 
FILE(GLOB SRCs SRCs/*.cxx)


## 2.1. VTK   
list(APPEND LIBS vtkPVCatalyst vtkPVVTKExtensionsDefault vtkPVClientServerCoreCore)
list(APPEND LIBS vtkFiltersParallelGeometry vtkFiltersParallelMPI vtkFiltersParallel)
list(APPEND LIBS vtkInteractionStyle vtkWrappingPythonCore)  
find_package(ParaView 4.1 REQUIRED COMPONENTS ${LIBS})
include("${PARAVIEW_USE_FILE}")

include(${VTK_USE_FILE})
list(APPEND LIBS ${VTK_LIBRARIES})  

find_package(PythonLibs 2.7 EXACT REQUIRED)
include_directories("${PYTHON_INCLUDE_PATH}")


## 3.1. LIBRARY   
set(VTK_MY_LIB_NAME vtkXYZ)
set(VTK_MY_PYTHON_LIB_NAME ${VTK_MY_LIB_NAME}Python)
add_library(${VTK_MY_LIB_NAME}Lib ${SRCs}) ## ?? 
target_link_libraries(${VTK_MY_LIB_NAME}Lib) 
 

## 4.1. WRAPPING  
include(vtkWrapPython)
VTK_WRAP_PYTHON3(${VTK_MY_PYTHON_LIB_NAME} PythonSRCs "${SRCs}") 
add_library(${VTK_MY_PYTHON_LIB_NAME}D ${PythonSRCs} ${SRCs})
target_link_libraries(${VTK_MY_PYTHON_LIB_NAME}D ${VTK_LIBRARIES} vtkWrappingPythonCore ${VTK_PYTHON_LIBRARIES})

## 4.2.  
add_library(${VTK_MY_PYTHON_LIB_NAME} MODULE ${VTK_MY_PYTHON_LIB_NAME}Init.cxx)

set(VTK_MODULES_USED vtkInteractionStyle)
set(VTK_PYTHOND_LIBS)
foreach(TMP_LIB ${VTK_MODULES_USED})
  set(VTK_PYTHOND_LIBS ${VTK_PYTHOND_LIBS} ${TMP_LIB}PythonD)
endforeach()
target_link_libraries(${VTK_MY_PYTHON_LIB_NAME} ${VTK_MY_PYTHON_LIB_NAME}D ${VTK_PYTHOND_LIBS})


set_target_properties(${VTK_MY_PYTHON_LIB_NAME} PROPERTIES PREFIX "")
set_target_properties(${VTK_MY_PYTHON_LIB_NAME} PROPERTIES SUFFIX ".so")

## SEE :
##     https://gitlab.kitware.com/vtk/vtk/blob/master/Wrapping/Python/README.md
##     https://github.com/patmarion/PythonWrappedPluginExample
##     https://github.com/paulmelis/vtk-wrapping-example 
