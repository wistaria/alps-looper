add_definitions(-DLOOP_CONFIG_HEADER="extras/localsus/localsus_config.h")

add_executable(localsus ${PROJECT_SOURCE_DIR}/loop.C ${PROJECT_SOURCE_DIR}/loop_model.C
  ${PROJECT_SOURCE_DIR}/path_integral.C)
target_link_libraries(localsus ${ALPS_LIBRARIES})

add_alps_test(localsus)
