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
#include "topology/RingTopology.hpp"
#include "utility/cli.hpp"

Utility::cliArguments a;
std::vector<Utility::Particle> particles;
MPI_Datatype mpiParticleType;

std::shared_ptr<Simulation> createContext() {
    // create topology
    std::shared_ptr<RingTopology> ringTopology = std::make_shared<RingTopology>();

    // domain decomposition
    std::shared_ptr<AtomDecomposition> atomDecomposition = std::make_shared<AtomDecomposition>();

    // create potential... TODO: pass correct parameters epsilon, sigma, nu
    std::shared_ptr<LennardJones> lj = std::make_shared<LennardJones>(1.0, 1.0);
    std::shared_ptr<AxilrodTeller> axilrodTeller = std::make_shared<AxilrodTeller>(1.0);

    // create algorithm
    std::shared_ptr<EAUTA> eauta = std::make_shared<EAUTA>();

    // set up simulation
    std::shared_ptr<Simulation> simulation =
        std::make_shared<Simulation>(a.iterations, a.respaStepSize, eauta, ringTopology, lj, axilrodTeller,
                                     atomDecomposition, &mpiParticleType, particles, a.deltaT, a.gForce, a.outputCSV);
    return simulation;
}

int main(int argc, char* argv[]) {
    // init MPI
    MPI_Init(&argc, &argv);

    int worldSize, worldRank;
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);

    // parse cli arguments
    std::vector<std::string> args;

    for (int i = 1; i < argc; i++) {
        args.push_back(argv[i]);
    }

    a = Utility::cliParse(args);

    // create particleMPIType
    mpiParticleType = Utility::Particle::GetMPIType();
    MPI_Type_commit(&mpiParticleType);

    // load particle input data
    Utility::getParticlesFromCSV(a.inputCSV, particles);

    std::shared_ptr<Simulation> simulation = createContext();

    simulation->Init();

    // execute simulation
    simulation->Start();

    // finalize
    MPI_Type_free(&mpiParticleType);
    MPI_Finalize();

    return 0;
}