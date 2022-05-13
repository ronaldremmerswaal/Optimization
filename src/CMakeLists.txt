add_library(optim
  m_linesearch.f90
  m_multid_optim.f90
  m_oned_rootfinding.f90
  m_poly_roots.f90
  m_optimization.f90
  )

set_target_properties (optim PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

target_include_directories(optim PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
