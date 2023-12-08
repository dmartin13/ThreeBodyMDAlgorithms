cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(
  3BMDA
  VERSION 1.0
  LANGUAGES C CXX
  DESCRIPTION "Three Body MD Algorithms")

# so we can see which compile flags are used by gcc
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(TESTS_3BMDA "Build test programs" OFF)
option(PROFILE_3BMDA "Profile mode" OFF)
option(USE_OMP "Use OPENMP" OFF)
option(MEASURESIMSTEP_3BMDA "Measure one timestep" OFF)
option(VLEVEL "Verbosity Level" -1)

find_package(Eigen3 REQUIRED)
find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RELEASE)
endif()

if(VLEVEL GREATER_EQUAL 0)
  add_compile_definitions(VLEVEL=${VLEVEL})
endif()

if(CMAKE_COMPILER_IS_GNUCXX)
  message(STATUS "GCC detected, adding compile flags")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
  if(${OpenMP_FOUND})
    message(STATUS "OpenMP detected, adding compiler definition")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
    add_compile_definitions(OPENMPAVAIL)
  else()
    unset(OPENMPAVAIL)
  endif(${OpenMP_FOUND})
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -pg -Og")
  set(CMAKE_CXX_FLAGS_RELEASE
      "${CMAKE_CXX_FLAGS_RELEASE} -march=native -fno-omit-frame-pointer")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO " -pg -Og -march=native")
endif(CMAKE_COMPILER_IS_GNUCXX)

add_subdirectory(src)
# add_subdirectory(tools)

# MAIN

include_directories(SYSTEM ${MPI_INCLUDE_PATH})

file(GLOB_RECURSE MAIN_SOURCES CONFIGURE_DEPENDS "src/*.hpp" "src/*.cpp")

# set(MAIN_SOURCES src/main.cpp src/topology/Topology.cpp
# src/topology/RingTopology.cpp src/topology/CartTopology.cpp
# src/decomposition/AtomDecomposition.cpp
# src/decomposition/RegularGridDecomposition.cpp
# src/decomposition/DomainDecomposition.cpp src/algorithm/Algorithm.cpp
# src/algorithm/NATA.cpp src/algorithm/AUTA.cpp src/algorithm/P3BCA.cpp
# src/utility/utility.cpp src/potential/AxilrodTeller.cpp
# src/potential/Potential.cpp src/potential/LennardJones.cpp
# src/simulation/Simulation.cpp src/MPIReporter.cpp)

add_executable(main main.cpp ${MAIN_SOURCES})
target_include_directories(main PRIVATE src)

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

# target_include_directories( main PUBLIC
# ${CMAKE_SOURCE_DIR}/src/external/rapidjson/include)

# TOOLS
file(GLOB_RECURSE TOOLS_SOURCES CONFIGURE_DEPENDS "tools/*.hpp" "tools/*.cpp")

add_executable(particlegenerator particleToCSV.cpp ${TOOLS_SOURCES})

target_include_directories(particlegenerator PRIVATE tools)

set_target_properties(
  particlegenerator PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                               ${CMAKE_BINARY_DIR}/executables/tools)

# TESTS

if(TESTS_3BMDA)

  include(GoogleTest)

  set(INSTALL_GTEST OFF)
  add_subdirectory(external/gtest)

  file(
    GLOB_RECURSE
    TEST_SOURCES
    CONFIGURE_DEPENDS
    "src/*.hpp"
    "src/*.cpp"
    "tests/*.hpp"
    "tests/*.cpp"
    "tools/*.hpp"
    "tools/*.cpp")

  # set(TEST_SOURCES tests/tests.cpp src/topology/Topology.cpp
  # src/topology/RingTopology.cpp src/topology/CartTopology.cpp
  # src/decomposition/AtomDecomposition.cpp
  # src/decomposition/RegularGridDecomposition.cpp
  # src/decomposition/DomainDecomposition.cpp src/algorithm/Algorithm.cpp
  # src/algorithm/NATA.cpp src/algorithm/AUTA.cpp src/algorithm/P3BCA.cpp
  # src/utility/utility.cpp src/potential/AxilrodTeller.cpp
  # src/potential/Potential.cpp src/simulation/Simulation.cpp
  # src/tools/ClosestPackedGenerator.cpp src/tools/ClusteredGaussGenerator.cpp
  # src/tools/GaussGenerator.cpp src/tools/GridGenerator.cpp
  # src/tools/ParticleGenerator.cpp src/tools/UniformGenerator.cpp
  # src/utility/utility.cpp)

  add_executable(tests ${TEST_SOURCES})

  target_include_directories(tests PRIVATE src)
  target_include_directories(tests PRIVATE external/gtest/include)

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
    tests PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                                ${CMAKE_BINARY_DIR}/executables/tests)

endif(TESTS_3BMDA)
