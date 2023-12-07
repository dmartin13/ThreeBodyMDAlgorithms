include_directories(SYSTEM ${MPI_INCLUDE_PATH})

add_subdirectory(utility)
add_subdirectory(tools)
add_subdirectory(external/rapidyaml)

set(MAIN_SOURCES
    main.cpp
    topology/Topology.cpp
    topology/RingTopology.cpp
    topology/CartTopology.cpp
    decomposition/AtomDecomposition.cpp
    decomposition/RegularGridDecomposition.cpp
    decomposition/DomainDecomposition.cpp
    algorithm/Algorithm.cpp
    algorithm/NATA.cpp
    algorithm/AUTA.cpp
    algorithm/P3BCA.cpp
    utility/utility.cpp
    potential/AxilrodTeller.cpp
    potential/Potential.cpp
    potential/LennardJones.cpp
    simulation/Simulation.cpp
    MPIReporter.cpp)

add_executable(main ${MAIN_SOURCES})
target_include_directories(main PRIVATE utility algorithm decomposition
                                        topology potential simulation)

target_link_libraries(main ${MPI_C_LIBRARIES})
target_link_libraries(main Eigen3::Eigen)
target_link_libraries(main OpenMP::OpenMP_CXX)

if(PROFILE_3BMDA)
  target_compile_definitions(main PUBLIC PROFILE_3BMDA)
endif()
if(USE_OMP)
  message(STATUS "add USE_OMP compile definition")
  target_compile_definitions(main PUBLIC USE_OMP)
endif()
if(MEASURESIMSTEP_3BMDA)
  target_compile_definitions(main PUBLIC MEASURESIMSTEP_3BMDA)
endif()

set_target_properties(main PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                                      ${CMAKE_BINARY_DIR}/executables)

target_include_directories(
  main PUBLIC ${CMAKE_SOURCE_DIR}/src/external/rapidjson/include)

# TESTS

if(TESTS_3BMDA)

  include(GoogleTest)

  set(INSTALL_GTEST OFF)
  add_subdirectory(external/gtest)
  include_directories(external/gtest/include)

  set(TEST_SOURCES
      tests/tests.cpp
      topology/Topology.cpp
      topology/RingTopology.cpp
      topology/CartTopology.cpp
      decomposition/AtomDecomposition.cpp
      decomposition/RegularGridDecomposition.cpp
      decomposition/DomainDecomposition.cpp
      algorithm/Algorithm.cpp
      algorithm/NATA.cpp
      algorithm/AUTA.cpp
      algorithm/P3BCA.cpp
      utility/utility.cpp
      potential/AxilrodTeller.cpp
      potential/Potential.cpp
      simulation/Simulation.cpp
      tools/ClosestPackedGenerator.cpp
      tools/ClusteredGaussGenerator.cpp
      tools/GaussGenerator.cpp
      tools/GridGenerator.cpp
      tools/ParticleGenerator.cpp
      tools/UniformGenerator.cpp
      utility/utility.cpp)

  add_executable(testrunner tests/testrunner.cpp)

  add_executable(tests ${TEST_SOURCES})
  target_link_libraries(tests gtest)
  target_link_libraries(tests ${MPI_C_LIBRARIES})
  target_link_libraries(tests Eigen3::Eigen)

  target_compile_definitions(tests PUBLIC TESTS_3BMDA)
  if(USE_OMP)
    message(STATUS "add USE_OMP compile definition")
    target_compile_definitions(tests PUBLIC USE_OMP)
  endif()

  # https://stackoverflow.com/questions/71971858/ctest-and-mpi-parallel-tests
  # solving the problem when gtest executes tests without mpi set_property(
  # TARGET tests PROPERTY CROSSCOMPILING_EMULATOR '${MPIEXEC_EXECUTABLE}
  # ${MPIEXEC_NUMPROC_FLAG} 1') gtest_discover_tests(tests TEST_PREFIX "tests/"
  # DISCOVERY_TIMEOUT 60)

  set_target_properties(
    tests testrunner PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                                ${CMAKE_BINARY_DIR}/executables/tests)

endif(TESTS_3BMDA)

# BENCHMARK
if(BENCHMARK_3BMDA)
  set(BENCHMARK_ENABLE_INSTALL OFF)
  set(BENCHMARK_ENABLE_TESTING OFF)
  set(BENCHMARK_INSTALL_DOCS OFF)
  set(BENCHMARK_DOWNLOAD_DEPENDENCIES ON)

  add_subdirectory(external/benchmark)

  add_executable(benchrunner benchmarks/benchrunner.cpp)

  set(BENCHMARK_SOURCES
      benchmarks/benchmain.cpp
      benchmarks/MPIBenchmark.cpp
      benchmarks/AUTAContext.cpp
      benchmarks/NATAContext.cpp
      benchmarks/P3BCAContext.cpp
      benchmarks/Context.cpp
      benchmarks/SingleIteration.cpp
      benchmarks/SingleIterationOnlySimStep.cpp
      topology/Topology.cpp
      topology/RingTopology.cpp
      topology/CartTopology.cpp
      decomposition/AtomDecomposition.cpp
      decomposition/RegularGridDecomposition.cpp
      decomposition/DomainDecomposition.cpp
      algorithm/Algorithm.cpp
      algorithm/NATA.cpp
      algorithm/AUTA.cpp
      algorithm/P3BCA.cpp
      utility/utility.cpp
      potential/AxilrodTeller.cpp
      potential/Potential.cpp
      simulation/Simulation.cpp
      tools/ClosestPackedGenerator.cpp
      tools/ClusteredGaussGenerator.cpp
      tools/GaussGenerator.cpp
      tools/GridGenerator.cpp
      tools/ParticleGenerator.cpp
      tools/UniformGenerator.cpp)

  add_executable(benchmain ${BENCHMARK_SOURCES})

  target_include_directories(benchmain PRIVATE thirdparty/rapidyaml/src)

  target_link_libraries(benchmain benchmark)
  target_link_libraries(benchmain ${MPI_C_LIBRARIES})
  target_link_libraries(benchmain Eigen3::Eigen)
  target_link_libraries(benchmain ryml)

  target_compile_definitions(benchmain PUBLIC BENCHMARK_3BMDA)
  if(USE_OMP)
    message(STATUS "add USE_OMP compile definition")
    target_compile_definitions(benchmain PUBLIC USE_OMP)
  endif()

  set(BENCHMARK_CLUSTER_SOURCES
      benchmarks/benchmain_cluster.cpp
      benchmarks/MPIBenchmark.cpp
      benchmarks/AUTAContext.cpp
      benchmarks/NATAContext.cpp
      benchmarks/P3BCAContext.cpp
      benchmarks/Context.cpp
      benchmarks/SingleIteration.cpp
      benchmarks/SingleIterationOnlySimStep.cpp
      topology/Topology.cpp
      topology/RingTopology.cpp
      topology/CartTopology.cpp
      decomposition/AtomDecomposition.cpp
      decomposition/RegularGridDecomposition.cpp
      decomposition/DomainDecomposition.cpp
      algorithm/Algorithm.cpp
      algorithm/NATA.cpp
      algorithm/AUTA.cpp
      algorithm/P3BCA.cpp
      utility/utility.cpp
      potential/AxilrodTeller.cpp
      potential/Potential.cpp
      simulation/Simulation.cpp
      tools/ClosestPackedGenerator.cpp
      tools/ClusteredGaussGenerator.cpp
      tools/GaussGenerator.cpp
      tools/GridGenerator.cpp
      tools/ParticleGenerator.cpp
      tools/UniformGenerator.cpp)

  add_executable(benchmain_cluster ${BENCHMARK_CLUSTER_SOURCES})

  target_include_directories(benchmain_cluster PRIVATE thirdparty/rapidyaml/src)

  target_link_libraries(benchmain_cluster benchmark)
  target_link_libraries(benchmain_cluster ${MPI_C_LIBRARIES})
  target_link_libraries(benchmain_cluster Eigen3::Eigen)
  target_link_libraries(benchmain_cluster ryml)

  target_compile_definitions(benchmain_cluster PUBLIC BENCHMARK_3BMDA)
  if(USE_OMP)
    message(STATUS "add USE_OMP compile definition")
    target_compile_definitions(benchmain_cluster PUBLIC USE_OMP)
  endif()

  set_target_properties(
    benchrunner benchmain benchmain_cluster
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY
               ${CMAKE_BINARY_DIR}/executables/benchmark)

endif(BENCHMARK_3BMDA)
