add_definitions(-DLOOP_CONFIG_HEADER="extras/xyz/xyz_config.h")

add_library(xyz_core ${PROJECT_SOURCE_DIR}/loop_model.C ${PROJECT_SOURCE_DIR}/diag.C
  ${PROJECT_SOURCE_DIR}/path_integral.C)
target_link_libraries(xyz_core ${ALPS_LIBRARIES})

add_executable(xyz ${PROJECT_SOURCE_DIR}/loop.C)
target_link_libraries(xyz xyz_core)

if(LOOPER_BUILD_TESTS)
  add_executable(xyz_ns ${PROJECT_SOURCE_DIR}/loop_ns.C)
  target_link_libraries(xyz_ns xyz_core)
  add_alps_test(xyz_ns)
  configure_file(xyz_model.xml.in xyz_model.xml)
endif(LOOPER_BUILD_TESTS)
