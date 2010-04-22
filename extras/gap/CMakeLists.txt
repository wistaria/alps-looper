add_definitions(-DLOOP_CONFIG_HEADER="extras/gap/gap_config.h")

add_library(gap_core ${PROJECT_SOURCE_DIR}/loop_model.C
  ${PROJECT_SOURCE_DIR}/path_integral.C)
target_link_libraries(gap_core ${ALPS_LIBRARIES})

add_executable(gap ${PROJECT_SOURCE_DIR}/loop.C)
target_link_libraries(gap gap_core)

if(LOOPER_BUILD_TESTS)
  add_executable(gap_ns ${PROJECT_SOURCE_DIR}/loop_ns.C)
  target_link_libraries(gap_ns gap_core)
  add_alps_test(gap_ns)
endif(LOOPER_BUILD_TESTS)
