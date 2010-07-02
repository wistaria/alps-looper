add_definitions(-DLOOP_CONFIG_HEADER="extras/transmag/transmag_config.h")

add_executable(transmag ${PROJECT_SOURCE_DIR}/loop.C ${PROJECT_SOURCE_DIR}/loop_model.C
  ${PROJECT_SOURCE_DIR}/path_integral.C)
target_link_libraries(transmag ${ALPS_LIBRARIES})

if(LOOPER_BUILD_TESTS)
  add_alps_test(transmag)
  configure_file(alternating_chain.xml.in alternating_chain.xml)
endif(LOOPER_BUILD_TESTS)
