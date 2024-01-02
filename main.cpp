#include <mpi.h>

#include <iostream>
#include <string>
#include <vector>

#include "MPIReporter.hpp"
#include "algorithm/EAUTA.hpp"
#include "decomposition/AtomDecomposition.hpp"
#include "potential/AxilrodTeller.hpp"
#include "potential/LennardJones.hpp"
#include "simulation/Simulation.hpp"
#include "tests/helpers/axilrodTellerMuto.hpp"
#include "tests/helpers/lennardJones.hpp"
#include "topology/RingTopology.hpp"
#include "utility/cli.hpp"

Utility::cliArguments a;
std::vector<Utility::Particle> particles;
MPI_Datatype mpiParticleType;

int main(int argc, char* argv[]) { return 0; }