cmake_minimum_required(VERSION 3.13)

project(optim C Fortran) 

# Unit testing
set(ENABLE_UNIT_TESTS_OPTIM False CACHE BOOL "build and run tests")
if (ENABLE_UNIT_TESTS_OPTIM)
  add_subdirectory(tests)
endif ()

# Add directiories
add_subdirectory(src)
