cmake_minimum_required(VERSION 3.5.1)

## mio, levmar, etc
include_directories(${MODULES_DIR})

## Setup Release and Debug variables
include("${MODULES_DIR}/mio/cmake/DefaultConfigTypes.cmake")

## exponential test
add_executable(exponential_test "${CMAKE_CURRENT_LIST_DIR}/exponential.cpp")
target_link_libraries(exponential_test levmar lapack blas)

## trig test
add_executable(trig_test "${CMAKE_CURRENT_LIST_DIR}/trig.cpp")
target_link_libraries(trig_test levmar lapack blas)