#include <mpi.h>

#include <format>
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

std::shared_ptr<Simulation> createContext() {
    // create topology
    std::shared_ptr<RingTopology> ringTopology = std::make_shared<RingTopology>();

    // domain decomposition
    std::shared_ptr<AtomDecomposition> atomDecomposition = std::make_shared<AtomDecomposition>();

    // create potential
    std::shared_ptr<LennardJones> lj = std::make_shared<LennardJones>(a.epsilon, a.sigma);
    std::shared_ptr<AxilrodTeller> axilrodTeller = std::make_shared<AxilrodTeller>(a.nu);

    // create algorithm
    std::shared_ptr<EAUTA> eauta = std::make_shared<EAUTA>();

    // set up simulation
    std::shared_ptr<Simulation> simulation =
        std::make_shared<Simulation>(a, a.iterations, a.respaStepSize, eauta, ringTopology, lj, axilrodTeller,
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

void forceCalcTest() {
    const double sigma = 1;
    const double epsilon = 1;
    const double nu = 1;

    const auto p0 = Eigen::Vector3d{1.0, 1.0, 0.0};
    const auto p1 = Eigen::Vector3d{5.0, 9.0, 0.0};
    const auto p2 = Eigen::Vector3d{9.0, 1.0, 0.0};

    std::array<Eigen::Vector3d, 3> pairwiseForces = {Eigen::Vector3d{0, 0, 0}, Eigen::Vector3d{0, 0, 0},
                                                     Eigen::Vector3d{0, 0, 0}};

    std::array<Eigen::Vector3d, 3> triwiseForces = {Eigen::Vector3d{0, 0, 0}, Eigen::Vector3d{0, 0, 0},
                                                    Eigen::Vector3d{0, 0, 0}};

    // calculate pairwise interactions
    auto fp0p1 = calculateLJForce(p0, p1, sigma, epsilon);
    auto fp0p2 = calculateLJForce(p0, p2, sigma, epsilon);
    auto fp1p2 = calculateLJForce(p1, p2, sigma, epsilon);

    pairwiseForces[0] += std::get<0>(fp0p1);
    pairwiseForces[0] += std::get<0>(fp0p2);

    pairwiseForces[1] += std::get<1>(fp0p1);
    pairwiseForces[1] += std::get<0>(fp1p2);

    pairwiseForces[2] += std::get<1>(fp0p2);
    pairwiseForces[2] += std::get<1>(fp1p2);

    // calculate triplet interactions
    auto fp0p1p2 = calculateATMForce(p0, p1, p2, nu);

    triwiseForces[0] += std::get<0>(fp0p1p2);

    triwiseForces[1] += std::get<1>(fp0p1p2);

    triwiseForces[2] += std::get<2>(fp0p1p2);

    for (int i = 0; i < 3; ++i) {
        std::cout << "total two-body force of p" << i << ": (" << pairwiseForces[i][0] << ", " << pairwiseForces[i][1]
                  << ", " << pairwiseForces[i][2] << "), total three-body force of p" << i << ": ("
                  << triwiseForces[i][0] << ", " << triwiseForces[i][1] << ", " << triwiseForces[i][2]
                  << "), total sum: (" << std::format("{}", pairwiseForces[i][0] + triwiseForces[i][0]) << ", "
                  << std::format("{}", pairwiseForces[i][1] + triwiseForces[i][1]) << ", "
                  << std::format("{}", pairwiseForces[i][2] + triwiseForces[i][2]) << ")" << std::endl;
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

    // if (worldRank == 0) {
    //     forceCalcTest();
    // }

    // finalize
    MPI_Type_free(&mpiParticleType);
    MPI_Finalize();

    return 0;
}