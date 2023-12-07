add_executable(
  particlegenerator
  particleToCSV.cpp
  ParticleGenerator.cpp
  ClusteredGaussGenerator.cpp
  ClosestPackedGenerator.cpp
  GaussGenerator.cpp
  GridGenerator.cpp
  UniformGenerator.cpp)

set_target_properties(
  particlegenerator PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                               ${CMAKE_BINARY_DIR}/executables/tools)
