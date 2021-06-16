set(TESTS alternating_tensor bond_matrix compatibility edge_type
    filtered_graph find_bridge generate_01 gesvd lattice llsp
    model_parameter omp_chunksize padding parse_vector permutation poisson_distribution
    postorder prime_factorization sign site_matrix spinmodel_helper syev
    temperature virtual_lattice weight)

foreach(name ${TESTS})
  add_executable(${name} ${name}.C)
  target_link_libraries(${name} alps Boost::mpi ${ALPS_EXTRA_LIBRARIES})
  add_alps_test(${name})
endforeach(name)

set(TESTS atomic union_find)

set(ATOMIC_SRC "")
if(CMAKE_ASM)
  enable_language(ASM)
  set_property(SOURCE looper/atomic.s PROPERTY LANGUAGE ASM)
  set(ATOMIC_SRC ${PROJECT_SOURCE_DIR}/looper/atomic.s)
endif(CMAKE_ASM)

foreach(name ${TESTS})
  add_executable(${name} ${name}.C ${ATOMIC_SRC})
  target_link_libraries(${name} alps Boost::mpi ${ALPS_EXTRA_LIBRARIES})
  add_alps_test(${name})
endforeach(name)
