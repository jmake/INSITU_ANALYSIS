## 
cmake_minimum_required(VERSION 3.0)

project(PLEPP VERSION 0.1 LANGUAGES C CXX Fortran)  
set(CMAKE_CXX_STANDARD 11)

set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Release")

set(PLEPP_ROOT "../../LIBPLEPP/" CACHE FILEPATH "path to src and include")

## Include Directories (for all targets)
include_directories(${PLEPP_ROOT}/include) #${CMAKE_BINARY_DIR})

## Wildcard all the sources in src
file(GLOB PLEPP_SOURCES ${PLEPP_ROOT}/src/*.cpp ${PLEPP_ROOT}/src/*.c)
add_library(PLEPP ${PLEPP_SOURCES})  

## Testing   
## 01 
set(exec01 iam_both.x)  
add_executable(${exec01} dummy.cxx)  
target_link_libraries(${exec01} PLEPP)  

## 03   
set(exec03 iam_alya.x)  
add_executable(${exec03} dummy.f90)  
set_property(TARGET ${exec03} PROPERTY LINKER_LANGUAGE Fortran)
target_link_libraries(${exec03} PLEPP)

## 02   
set(exec02 iam_hemelb.x)
add_executable(${exec02} dummy.cxx)
target_compile_definitions(${exec02} PRIVATE -DWITHALYA)
target_link_libraries(${exec02} PLEPP)


enable_testing() ## SEE [1]

## EXAMPLE01   
set(case01  "mpmd1")
#set(torun01 "mpirun -np 2 ./${exec01} HEMELB : -np 1 ./${exec01} ALYA")
##set(fname01 "${case01}.sh")
##file(WRITE ${fname01} ${torun01})
#add_test(NAME ${case01} COMMAND mpirun "-np 2" "./${exec01}" "HEMELB" ":" "-np 1" "./${exec01}" "ALYA")
#add_test(NAME ${case01} COMMAND mpirun -np 2 ./${exec01} HEMELB : -np 1 ./${exec01} ALYA)


## EXAMPLE02 
#set(case01  "mpmd2")
#set(fname01 "${case01}.sh")
#set(torun01 " 0-1 ./${exec01} HEMELB \n 2-3 ./${exec01} ALYA \n 4-7 ./${exec01} ALYA")
#file(WRITE ${fname01} ${torun01})
#add_test(NAME ${case01} COMMAND "srun" "--multi-prog" ${fname01})


## EXAMPLE02 
#set(case01  "mpmd3")
#set(fname01 "${case01}.sh")
#set(torun01 " 0-1 ./${exec01} HEMELB \n 2-3 ./${exec03} ALYA \n 4-7 ./${exec03} ALYA")
#file(WRITE ${fname01} ${torun01})
#add_test(NAME ${case01} COMMAND "srun" "--multi-prog" ${fname01})


##
## [1] vim Testing/Temporary/LastTest.log  
## 
