add_definitions(-DLOOP_CONFIG_HEADER="extras/gap/gap_config.h")

add_executable(gap ${PROJECT_SOURCE_DIR}/loop.C ${PROJECT_SOURCE_DIR}/loop_model.C
  ${PROJECT_SOURCE_DIR}/path_integral.C)
target_link_libraries(gap ${ALPS_LIBRARIES})

if(LOOPER_BUILD_TESTS)
  add_alps_test(gap)
endif(LOOPER_BUILD_TESTS)
