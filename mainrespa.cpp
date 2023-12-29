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

void gatherAndPrintMessages() {
    int worldSize, worldRank;
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);

    std::vector<std::string> allMyMessages = MPIReporter::instance()->GetAllMessages();
    int numMyMessages = allMyMessages.size();

    std::vector<int> numAllMessages;
    numAllMessages.resize(worldSize);

    MPI_Gather(&numMyMessages, 1, MPI_INT, numAllMessages.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (worldRank == 0) {
        std::vector<std::string> allMessages;

        for (std::string& str : allMyMessages) {
            allMessages.push_back(str);
        }

        for (int i = 1; i < worldSize; i++) {
            for (int j = 0; j < numAllMessages[i]; j++) {
                MPI_Status status;
                int numRecv;

                MPI_Probe(i, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                MPI_Get_count(&status, MPI_CHAR, &numRecv);

                char* buf = new char[numRecv];

                MPI_Recv(buf, numRecv, MPI_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                std::string str(buf, numRecv);

                allMessages.push_back(str);

                delete[] buf;
            }
        }

        for (std::string& m : allMessages) {
            std::cout << m << std::endl;
        }

    } else {
        std::vector<MPI_Request> requests;
        requests.resize(numMyMessages);
        for (int j = 0; j < numMyMessages; j++) {
            // MPI_Request req;
            // requests[j] = req;

            MPI_Isend(allMyMessages[j].c_str(), allMyMessages[j].length(), MPI_CHAR, 0, 0, MPI_COMM_WORLD,
                      &(requests[j]));
        }
        MPI_Waitall(numMyMessages, requests.data(), MPI_STATUSES_IGNORE);
    }
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

    // print messages
    gatherAndPrintMessages();

    // finalize
    MPI_Type_free(&mpiParticleType);
    MPI_Finalize();

    return 0;
}