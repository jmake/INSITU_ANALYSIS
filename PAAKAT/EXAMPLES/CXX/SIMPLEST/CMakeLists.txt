cmake_minimum_required(VERSION 3.0) 
project(PLEPP VERSION 0.1 LANGUAGES C CXX)  
set(CMAKE_CXX_STANDARD 11)


## VTK   
list(APPEND LIBS vtkPVCatalyst vtkPVVTKExtensionsDefault vtkPVClientServerCoreCore)
list(APPEND LIBS vtkFiltersParallelGeometry vtkFiltersParallelMPI vtkFiltersParallel)
#list(APPEND LIBS vtkPVPythonCatalyst) # Only for python 
#set(CMAKE_CXX_FLAGS "-DWITH_PYTHON")  # Only for python 
find_package(ParaView 4.1 REQUIRED COMPONENTS ${LIBS})
include("${PARAVIEW_USE_FILE}")


## Testing   
## 01 
set(exec01 test_cxx.x)  
add_executable(${exec01} test.cxx)  
set_property(TARGET ${exec01} PROPERTY LINKER_LANGUAGE CXX)
target_link_libraries(${exec01} ${LIBS})   

