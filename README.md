# Overview

# Cloning the project
Please use the git parameter ```--recurse-submodules``` when cloning the project, so all submodules are loaded.

# Dependencies
- [cmake](https://cmake.org/) >= 3.10 
- [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) 3.4
- [OpenMPI](https://www.open-mpi.org/), [MPICH](https://www.mpich.org/) or Intel MPI

This program was tested with cmake 3.16.3, eigen 3.4.0, OpenMPI 4.1.2 and gcc 13.1

# Building the Project
Execute the following command in the root directory of the project: ```mkdir build && cd build && cmake ..```. This prepares the cmake project for compilation.

## Simple build without additional options
Execute `make` in the build directory after executing the above command. The main program is located in build/executables. This command will also build the particle generator into build/executables/tools.

# Using the Program

## Parameters for the main program
- `-a`: The algorithm to be used for the simulation
  - Possible values: `NATA`, `AUTA`, `P3BCA`
- `-i`: Number of simulation steps (iterations)
- `-d`: Delta-time (duration of one time-step)
- `-gx`: Gravity in x-direction (analogous: `-gy`, `-gz`)
- `-c`: The cutoff distance in the case of the p3bca algorithm
- `-csv`: relative path to a csv-file containing particles, generated with the particle generator
- `-decomp`: The decomposition strategy for the p3bca algorithm.
  - possible values: `optimal`, `naive`
  - `optimal`: A dynamic 1D, 2D or 3D regular grid decomposition based on the prime decomposition of the number of processors
  - `naive`: Subdivides the domain along the x-axis in equal slices
- `-o`: specifies a csv outputfile, that stores all updated particles after the simulation

## Example Execution

## Simple measurement of one Simulation Step
- Pass the following parameter to CMake: ```cmake .. -DMEASURESIMSTEP_3BMDA=ON```
- The main program can be executed as in the upper step. As output, the times of all processors, as well as the maximum time over all processors are printed on the console in a JSON format (in nanoseconds).

## Particle Generator
The following particle distributions can be generated:
- Uniform: Generates uniformly distributed particles
  - Parameters
    - `--numparticles`: Number of particles to be generated
    - `--blx`: Domain-size in x-dir (analogous: `--bly`, `--blz`)
  - Example: ```./particlegenerator --generator uniform --numparticles 2048 --output uniform.csv --vx 0 --vy 0 --vz 0 --mass 1 --seed0 926762934 --blx 10 --bly 10 --blz 10```
- Grid: Generates particles with equal spacing
  - Parameters
    - `--particlesperx`: number of particles in x-dir (analogous: `--particlespery`, `--particlesperz`)
    - `--particlespacing`: spacing between particles in all directions
  - Example: ```./particlegenerator --generator grid --output grid.csv --vx 0 --vy 0 --vz 0 --mass 1 --particlesperx 10 --particlespery 10 --particlesperz 10 --particlespacing 1```
- Closest Packed: Generates hexagonally closest packed particles
  - Parameters
    - `--blx`: Domain-size in x-dir (analogous: `--bly`, `--blz`)
    - `--particlespacing`: spacing between particles in all directions
  - Example: ```./particlegenerator --generator closestpacked --output cpacked.csv --vx 0 --vy 0 --vz 0 --mass 1 --blx 10 --bly 10 --blz 10 --particlespacing 0.5```
- Clustered Gauss: Generates uniform distributed gaussian clouds
  - Parameters
    - `--numparticles`: Total number of particles to be generated
    - `--blx`: Domain-size in x-dir (analogous: `--bly`, `--blz`)
    - `--numclusters`: Number of gaussian clouds to be generated
    - `--distmean`: Mean for the normal distribution in all directions of all gaussian clouds
    - `--diststddev`: Variance for the normal distribution in all directions of all gaussian clouds
  - Example: ```./particlegenerator --generator clusteredgauss --numparticles 2500 --distmean 0 --diststddev 1 --output cgauss.csv --vx 0 --vy 0 --vz 0 --mass 1 --seed0 926762934 --seed1 89347587 --blx 10 --bly 10 --blz 10 --numclusters 25```
- Gauss: Generates one gaussian cloud
  - Parameters
    - `--numparticles`: Total number of particles to be generated
    - `--blx`: Domain-size in x-dir (analogous: `--bly`, `--blz`)
    - `--distmean`: Mean for the normal distribution in all directions
    - `--diststddev`: Variance for the normal distribution  in all directions
  - Example: ```./particlegenerator --generator gauss --numparticles 2500 --distmean 5 --diststddev 1 --output gauss.csv --vx 0 --vy 0 --vz 0 --mass 1 --seed0 926762934 --blx 10 --bly 10 --blz 10```

### Parameters for all distributions:
- `--generator`: The particle generator to be used
  - Possible values: `uniform`, `grid`, `closestpacked`, `clusteredgauss`, `gauss`
- `--output`: Output csv-file that stores the generated particles
- `--vx`: Initial particle velocity in x-dir
- `--vy`: Initial particle velocity in y-dir
- `--vz` : Initial particle velocity in z-dir
- `--mass`: particle mass
- `--seed0`: seed for the uniform, clusteredgauss and gauss distribution
- `--seed1`: additional seed for the clusteredgauss distribution

## Additional CMake Parameters
The following parameters are available for cmake:
- TESTS_3BMDA: This will compile the unit tests into build/executables/tests.
  - Usage: ```cmake .. -DTESTS_3BMDA=ON && make```
  - Run tests: go to build/executables/tests and execute ./testrunner
- PROFILE_3BMDA: Compiles the main program containing the compiler directives so that individual methods can be measured in more detail. A JSON file is generated as output, which stores the times across all processors spent in different methods.
  - Usage: ```cmake .. -DPROFILE_3BMDA=ON && make```
  - An exemplary execution looks as follows: `mpiexec -n 4 ./main -a auta -i 1 -d 0.001 -gx 0 -gy 0 -gz 0 -c 0 -csv uniform_2048.csv -decomp optimal -op "profile_auta_uniform_2048_4.json"`
    In the example, individual methods of the All Unique Triplets algorithm are measured and stored in the file `profile_auta_uniform_2048_4.json`. Within the output file the measured times are represented in different ways:
      - `avg time all proc`: The average time of this method over all substeps and all processors
      - `total time all proc`: The sum of the times of this method over all substeps and processors
      - `total avg time union proc`: The sum of the average times over all processors of each substep.
      - `total max time union proc`: The sum of the maximum values over all processors of each substep
      - `times per proc`: prints for each processor `i in {0 ... p-1}` the times of all substeps of the respective method as a list
- MEASURESIMSTEP_3BMDA: Compiles the main program so that the time of a simulation step can be measured. 
  - Usage: ```cmake .. -DMEASURESIMSTEP_3BMDA=ON && make```
  - Measure one simulation-step: see the above bullet point "Simple measurement of one Simulation Step"
- VLEVEL: More verbose output is printed to the console when a simulation step is executed
  - Usage: ```cmake .. -DVLEVEL=1 && make```, then the program can be used as described in "Example Execution".

All default values are either `OFF`, or `-1` if a parameter is not specified.

# External Tools
The following external tools are used for this project:
- [Google Test](https://github.com/google/googletest)
- [RapidCSV](https://github.com/d99kris/rapidcsv)
- [RapidJSON](https://rapidjson.org/)
- [RapidYAML](https://github.com/biojppm/rapidyaml)

# References
- [1]: P. Koanantakool and K. Yelick, "A Computation- and Communication-Optimal Parallel Direct 3-Body Algorithm," SC '14: Proceedings of the International Conference for High Performance Computing, Networking, Storage and Analysis, 2014, pp. 363-374, doi: 10.1109/SC.2014.35.
- [2]: Marcelli, Gianluca. The role of three-body interactions on the equilibrium and non-equilibrium properties of fluids from molecular simulation. University of Kent (United Kingdom), 2001.
- [3]: Martin, David. A Comparison of Three-body Algorithms for Molecular Dynamics Simulations. Technical University of Munich, Nov 2022