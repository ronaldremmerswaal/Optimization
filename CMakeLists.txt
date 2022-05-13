

cmake_minimum_required(VERSION 3.13)

project(optim Fortran)

# Unit testing
if (ENABLE_UNIT_TESTS)
  enable_testing()
endif ()


add_library(optim)

# Alias (Daniel Pfeiffer, Effective CMake) - this allows other
# projects that use Pkg as a subproject to find_package(Nmspc::Pkg)
# which does nothing because Pkg is already part of the project

add_library(optim::optim ALIAS optim)
set(OPTIM_LIBRARIES optim::optim CACHE STRING "OPTIM library target")


add_subdirectory(src)


# Wonton targets

install(TARGETS optim
  EXPORT OPTIM_LIBRARIES
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION include
  INCLUDES DESTINATION include
  )

#-----------------------------------------------------------------------------
# Prepare output for configuration files to be used by projects importing Wonton
#-----------------------------------------------------------------------------

# export targets

install(EXPORT OPTIM_LIBRARIES
  FILE OPTIMTargets.cmake
  NAMESPACE optim::
  EXPORT_LINK_INTERFACE_LIBRARIES
  DESTINATION lib/cmake/optim)
