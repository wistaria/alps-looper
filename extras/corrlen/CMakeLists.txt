add_definitions(-DLOOP_CONFIG_HEADER="extras/corrlen/corrlen_config.h")

add_executable(corrlen ${PROJECT_SOURCE_DIR}/loop.C ${PROJECT_SOURCE_DIR}/loop_model.C
  ${PROJECT_SOURCE_DIR}/path_integral.C)
target_link_libraries(corrlen ${ALPS_LIBRARIES})

if(LOOPER_BUILD_TESTS)
  add_alps_test(corrlen)
endif(LOOPER_BUILD_TESTS)
