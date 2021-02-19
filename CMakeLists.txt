cmake_minimum_required (VERSION 3.0)
project (Heavy-ion_EM_field)

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_HOME_DIRECTORY}")

set (CMAKE_CXX_FLAGS "-O3 -fopenmp -Wall")

add_subdirectory (src)
