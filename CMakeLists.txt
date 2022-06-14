cmake_minimum_required(VERSION 3.13)

project(optim C Fortran) 

# Unit testing
set(ENABLE_UNIT_TESTS_OPTIM False CACHE BOOL "build and run tests")
if (ENABLE_UNIT_TESTS_OPTIM)
  add_subdirectory(tests)
endif ()

# Default to Release build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
  message("-- Build type is set to ${CMAKE_BUILD_TYPE}")
endif()

message("-- Using ${CMAKE_Fortran_COMPILER_ID} Fortran compiler")

if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
  set(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -Wall")
  set(CMAKE_Fortran_FLAGS_DEBUG   "-g -traceback")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ip -xHOST")
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
  set(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -Wall")
  set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g3")
  set(CMAKE_Fortran_FLAGS_RELEASE "-Ofast -march=native")
endif()



# Add directiories
add_subdirectory(src)