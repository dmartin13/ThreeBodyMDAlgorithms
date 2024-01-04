cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(
  3BMDA
  VERSION 1.0
  LANGUAGES C CXX
  DESCRIPTION "Multiple Timestepping for Two and Three Body Forces")

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

# Particle Generator
file(GLOB_RECURSE TOOLS_SOURCES CONFIGURE_DEPENDS "tools/particlegenerator/*.hpp" "tools/particlegenerator/*.cpp")

add_executable(particlegenerator tools/particlegenerator/particlegenerator.cpp ${TOOLS_SOURCES})

target_include_directories(particlegenerator PRIVATE tools/particlegenerator)

set_target_properties(
  particlegenerator PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                               ${CMAKE_BINARY_DIR}/executables/tools)

# MAINRESPA

file(GLOB_RECURSE MAINRESPA_SOURCES CONFIGURE_DEPENDS "src/*.hpp" "src/*.cpp")

add_executable(mainrespa mainrespa.cpp ${MAINRESPA_SOURCES})

target_include_directories(mainrespa PRIVATE src)
target_include_directories(mainrespa PRIVATE ${MPI_INCLUDE_PATH})

target_link_libraries(mainrespa ${MPI_C_LIBRARIES})
target_link_libraries(mainrespa Eigen3::Eigen)
target_link_libraries(mainrespa OpenMP::OpenMP_CXX)

if(PROFILE_3BMDA)
  message(STATUS "add PROFILE compile definition")
  target_compile_definitions(mainrespa PUBLIC PROFILE_3BMDA)
endif()
if(USE_OMP)
  message(STATUS "add USE_OMP compile definition")
  target_compile_definitions(mainrespa PUBLIC USE_OMP)
endif()
if(MEASURESIMSTEP_3BMDA)
  message(STATUS "add MEASURESIMSTEP compile definition")
  target_compile_definitions(mainrespa PUBLIC MEASURESIMSTEP_3BMDA)
endif()

set_target_properties(mainrespa PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                                           ${CMAKE_BINARY_DIR}/executables)

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

  add_executable(tests ${TEST_SOURCES})

  target_include_directories(tests PRIVATE src)
  target_include_directories(tests PRIVATE external/gtest/include)
  target_include_directories(tests PRIVATE ${MPI_INCLUDE_PATH})

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

  set_target_properties(tests PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                                         ${CMAKE_BINARY_DIR}/executables/tests)

endif(TESTS_3BMDA)
