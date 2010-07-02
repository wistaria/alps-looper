add_definitions(-DLOOP_CONFIG_HEADER="extras/xyz/xyz_config.h")

add_executable(xyz ${PROJECT_SOURCE_DIR}/loop.C ${PROJECT_SOURCE_DIR}/loop_model.C
  ${PROJECT_SOURCE_DIR}/diag.C ${PROJECT_SOURCE_DIR}/path_integral.C)
target_link_libraries(xyz ${ALPS_LIBRARIES})

if(LOOPER_BUILD_TESTS)
  add_alps_test(xyz)
  configure_file(xyz_model.xml.in xyz_model.xml)
endif(LOOPER_BUILD_TESTS)
