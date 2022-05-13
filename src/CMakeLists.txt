# set the sources from which to compile the target library

set(optim_SOURCES
  m_linesearch.f90
  m_multid_optim.f90
  m_oned_rootfinding.f90
  m_optimization.f90
  )

target_sources(optim PRIVATE ${optim_SOURCES})

# What headers should it use and install

set(optim_HEADERS
  )

list(TRANSFORM optim_HEADERS PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/")

set_target_properties(optim PROPERTIES
  PUBLIC_HEADER "${optim_HEADERS}")

target_include_directories(optim PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${optim_BINARY_DIR}>
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
  $<INSTALL_INTERFACE:include>)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(dialect "-ffree-form -std=f2008 -fimplicit-none")
    set(bounds "-fbounds-check")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(dialect "-stand f08 -free -implicitnone")
    set(bounds "-check bounds")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
    set(dialect "-Mfreeform -Mdclchk -Mstandard -Mallocatable=03")
    set(bounds "-C")
endif()

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")
