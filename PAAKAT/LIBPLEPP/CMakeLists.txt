## 
cmake_minimum_required(VERSION 3.0)

project(PLEPP VERSION 0.1 LANGUAGES C CXX Fortran)  
set(CMAKE_CXX_STANDARD 11)

set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Release")

## Include Directories (for all targets)
include_directories(include ${CMAKE_BINARY_DIR})

## Wildcard all the sources in src
file(GLOB SOURCES src/*.cpp src/*.c)
add_library(COMMDOM ${SOURCES})

## Testing   
set(exec01 testc01_1.x)
add_executable(${exec01} wrappers/Cpp/explorer.cxx)   
target_link_libraries(${exec01} COMMDOM)  

#set(exec02 testf01_1.x)
#add_executable(${exec02} wrappers/Fortran/test01_1.f90)
#target_link_libraries(${exec02} COMMDOM)

set(exec03 testf02_1.x)
add_executable(${exec03} wrappers/Fortran/test02_1.f90)
target_link_libraries(${exec03} COMMDOM)

enable_testing() ## SEE [1]
set(test_parameters mpirun -np 2 ${exec01} ./HEMELB : -np 1 ${exec01} ./ALYA)
add_test(coupling21 ${test_parameters})

##
## [1] vim Testing/Temporary/LastTest.log  
## 
