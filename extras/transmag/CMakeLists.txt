add_definitions(-DLOOP_CONFIG_HEADER="extras/transmag/transmag_config.h")

add_library(transmag_core ${PROJECT_SOURCE_DIR}/loop_model.C
  ${PROJECT_SOURCE_DIR}/path_integral.C)
target_link_libraries(transmag_core ${ALPS_LIBRARIES})

add_executable(transmag ${PROJECT_SOURCE_DIR}/loop.C)
target_link_libraries(transmag transmag_core)

if(LOOPER_BUILD_TESTS)
  add_executable(transmag_ns ${PROJECT_SOURCE_DIR}/loop_ns.C)
  target_link_libraries(transmag_ns transmag_core)
  add_alps_test(transmag_ns)
  configure_file(alternating_chain.xml.in alternating_chain.xml)
endif(LOOPER_BUILD_TESTS)
