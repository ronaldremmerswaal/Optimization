find_package(PFUNIT REQUIRED)
enable_testing()

# tests
add_pfunit_ctest (optim_tests
  TEST_SOURCES oned.pf multid.pf
  LINK_LIBRARIES optim
  )
