#include <mpi.h>

#include <iostream>
#include <string>
#include <vector>

#include "MPIReporter.hpp"
#include "algorithm/AUTA.hpp"
#include "algorithm/NATA.hpp"
#include "algorithm/P3BCA.hpp"
#include "decomposition/AtomDecomposition.hpp"
#include "decomposition/RegularGridDecomposition.hpp"
#include "potential/AxilrodTeller.hpp"
#include "simulation/Simulation.hpp"
#include "topology/CartTopology.hpp"
#include "topology/RingTopology.hpp"
#include "utility/cli.hpp"
#include "utility/decompositions.hpp"

#ifdef PROFILE_3BMDA
#include <numeric>

#include "rapidjson/document.h"
#include "rapidjson/stringbuffer.h"
#include "rapidjson/writer.h"
#endif

Utility::cliArguments a;
std::vector<Utility::Particle> particles;
MPI_Datatype mpiParticleType;

std::shared_ptr<Simulation> createNATAContext(std::string csvOut)
{
    // create topology
    std::shared_ptr<RingTopology> ringTopology = std::make_shared<RingTopology>();

    // domain decomposition
    std::shared_ptr<AtomDecomposition> atomDecomposition = std::make_shared<AtomDecomposition>();

    // create potential
    std::shared_ptr<AxilrodTeller> axilrodTeller = std::make_shared<AxilrodTeller>(1.0);

    // create algorithm
    std::shared_ptr<NATA> nata = std::make_shared<NATA>();

    // set up simulation
    // int iterations, std::shared_ptr<Algorithm> algorithm, std::shared_ptr<Topology> topology,
    // std::shared_ptr<Potential> potential, std::shared_ptr<DomainDecomposition> decomposition,
    // MPI_Datatype* mpiParticleType, std::vector<Utility::Particle>& particles, double dt,
    // Eigen::Vector3d gForce
    std::shared_ptr<Simulation> simulation =
        std::make_shared<Simulation>(a.iterations, nata, ringTopology, axilrodTeller, atomDecomposition,
                                     &mpiParticleType, particles, a.deltaT, a.gForce, csvOut);
    return simulation;
}

std::shared_ptr<Simulation> createP3BCAContext(std::string csvOut, std::vector<int> decomposition)
{
    // create topology
    std::shared_ptr<CartTopology> cartTopology = std::make_shared<CartTopology>(decomposition);

    // domain decomposition
    std::shared_ptr<RegularGridDecomposition> regularGridDecomposition = std::make_shared<RegularGridDecomposition>();

    // create potential
    std::shared_ptr<AxilrodTeller> axilrodTeller = std::make_shared<AxilrodTeller>(1.0);

    // create algorithm
    std::shared_ptr<P3BCA> p3bca = std::make_shared<P3BCA>(a.cutoff);

    // set up simulation
    // int iterations, std::shared_ptr<Algorithm> algorithm, std::shared_ptr<Topology> topology,
    // std::shared_ptr<Potential> potential, std::shared_ptr<DomainDecomposition> decomposition,
    // MPI_Datatype* mpiParticleType, std::vector<Utility::Particle>& particles, double dt,
    // Eigen::Vector3d gForce
    std::shared_ptr<Simulation> simulation =
        std::make_shared<Simulation>(a.iterations, p3bca, cartTopology, axilrodTeller, regularGridDecomposition,
                                     &mpiParticleType, particles, a.deltaT, a.gForce, csvOut);
    return simulation;
}

std::shared_ptr<Simulation> createAUTAContext(std::string csvOut)
{
    // create topology
    std::shared_ptr<RingTopology> ringTopology = std::make_shared<RingTopology>();

    // domain decomposition
    std::shared_ptr<AtomDecomposition> atomDecomposition = std::make_shared<AtomDecomposition>();

    // create potential
    std::shared_ptr<AxilrodTeller> axilrodTeller = std::make_shared<AxilrodTeller>(1.0);

    // create algorithm
    std::shared_ptr<AUTA> auta = std::make_shared<AUTA>();

    // set up simulation
    // int iterations, std::shared_ptr<Algorithm> algorithm, std::shared_ptr<Topology> topology,
    // std::shared_ptr<Potential> potential, std::shared_ptr<DomainDecomposition> decomposition,
    // MPI_Datatype* mpiParticleType, std::vector<Utility::Particle>& particles, double dt,
    // Eigen::Vector3d gForce
    std::shared_ptr<Simulation> simulation =
        std::make_shared<Simulation>(a.iterations, auta, ringTopology, axilrodTeller, atomDecomposition,
                                     &mpiParticleType, particles, a.deltaT, a.gForce, csvOut);
    return simulation;
}

#ifdef PROFILE_3BMDA

std::string charToTimeUnit(const char& c)
{
    if (c == 0) {
        return std::string("nanosecond");
    } else if (c == 1) {
        return std::string("microsecond");
    } else if (c == 2) {
        return std::string("millisecond");
    } else if (c == 3) {
        return std::string("second");
    } else if (c == 4) {
        return std::string("minute");
    } else if (c == 5) {
        return std::string("hour");
    }
    return "";
}

void doTimingStuff(std::shared_ptr<Simulation> simulation, std::string outFile)
{
    // fetch profiling data
    std::map<std::string, std::pair<char, std::vector<int64_t>>> times = simulation->GetAlgorithm()->GetTimes();
    std::vector<double> hitrates = simulation->GetAlgorithm()->GetHitrates();
    // std::map<std::string, std::pair<char, double>> calcTimes = simulation->GetPotential()->GetAvgCalcTime();

    std::vector<float> valuesForSTVR;
    std::vector<CartRank> cartRanks;

    // init json object
    rapidjson::Document d;
    rapidjson::Value o(rapidjson::kObjectType);
    rapidjson::Value hr(rapidjson::kObjectType);
    rapidjson::Value hrPerProcRJ(rapidjson::kObjectType);
    rapidjson::Value stvrRJ(rapidjson::kArrayType);
    if (simulation->GetTopology()->GetWorldRank() == 0) {
        d.SetObject();
    }

    // gather all ranks
    {
        CartRank rank;
        std::shared_ptr<CartTopology> cartTopology = std::dynamic_pointer_cast<CartTopology>(simulation->GetTopology());
        if (cartTopology) {
            rank = cartTopology->GetCartRank();
        } else {
            rank = CartRank(simulation->GetTopology()->GetWorldRank());
        }
        std::array<int, 3> rankArray = rank.GetRank();

        std::vector<int> ranks;
        if (simulation->GetTopology()->GetWorldRank() == 0) {
            ranks.resize(simulation->GetTopology()->GetWorldSize() * rank.GetDimensions());
        }

        MPI_Gather(rankArray.data(), rank.GetDimensions(), MPI_INT, ranks.data(), rank.GetDimensions(), MPI_INT, 0,
                   simulation->GetTopology()->GetComm());

        if (simulation->GetTopology()->GetWorldRank() == 0) {
            for (int i = 0; i < simulation->GetTopology()->GetWorldSize(); i++) {
                if (rank.GetDimensions() == 1) {
                    cartRanks.push_back(CartRank(ranks[i * 1]));
                } else if (rank.GetDimensions() == 2) {
                    cartRanks.push_back(CartRank(ranks[i * 2], ranks[i * 2 + 1]));
                } else {
                    cartRanks.push_back(CartRank(ranks[i * 3], ranks[i * 3 + 1], ranks[i * 3 + 2]));
                }
            }
        }
    }

    // profile times calculation
    {
        for (const auto& [k, v] : times) {
            std::vector<int64_t> allElements;
            int numMyElements = v.second.size();
            std::vector<int> numAllElements;
            numAllElements.resize(simulation->GetTopology()->GetWorldSize());
            // elements from each process are gathered in order of their rank
            MPI_Gather(&numMyElements, 1, MPI_INT, numAllElements.data(), 1, MPI_INT, 0,
                       simulation->GetTopology()->GetComm());

            int maxNumElements = *(std::max_element(numAllElements.begin(), numAllElements.end()));

            std::vector<int> displacements;
            int sumDispl = 0;
            for (int i = 0; i < simulation->GetTopology()->GetWorldSize(); i++) {
                displacements.push_back(sumDispl);
                // check for over & underflow. https://stackoverflow.com/a/1514309
                if (maxNumElements > 0 && sumDispl > std::numeric_limits<int>::max() - maxNumElements) {
                    std::cout << "Overflow Warning for profiling in sumDispl" << std::endl;
                }
                if (maxNumElements < 0 && sumDispl < std::numeric_limits<int>::min() - maxNumElements) {
                    std::cout << "Underflow Warning for profiling in sumDispl" << std::endl;
                }
                sumDispl += maxNumElements;
            }

            allElements.resize(sumDispl);
            std::fill(allElements.begin(), allElements.end(), 0);

            MPI_Gatherv(v.second.data(), v.second.size(), MPI_INT64_T, allElements.data(), numAllElements.data(),
                        displacements.data(), MPI_INT64_T, 0, simulation->GetTopology()->GetComm());

            if (simulation->GetTopology()->GetWorldRank() == 0) {
                std::vector<int64_t> sumPerProc;
                sumPerProc.resize(simulation->GetTopology()->GetWorldSize());

                // rapidjson::Value timesPerProc(rapidjson::kArrayType);
                rapidjson::Value timesPerProc(rapidjson::kObjectType);

                for (int i = 0; i < simulation->GetTopology()->GetWorldSize(); i++) {
                    rapidjson::Value timesOfThisProc(rapidjson::kArrayType);
                    int64_t acc = 0;
                    for (int j = 0; j < maxNumElements; j++) {
                        if (j < numAllElements[i]) {
                            int64_t t = allElements[displacements[i] + j];

                            timesOfThisProc.PushBack(rapidjson::Value(t), d.GetAllocator());

                            // check for over & underflow. https://stackoverflow.com/a/1514309
                            if (t > 0 && acc > std::numeric_limits<int64_t>::max() - t) {
                                std::cout << "Overflow Warning for profiling" << std::endl;
                            }
                            if (t < 0 && acc < std::numeric_limits<int64_t>::min() - t) {
                                std::cout << "Underflow Warning for profiling" << std::endl;
                            }
                            acc += t;
#if defined(VLEVEL) && !defined(BENCHMARK_3BMDA) && !defined(TESTS_3BMDA) && VLEVEL > 0
                            std::string message = "proc " + std::to_string(i) + " added " + std::to_string(t) +
                                                  " measured timesteps for " + k;
                            MPIReporter::instance()->StoreMessage(simulation->GetTopology()->GetWorldRank(), message);
#endif
                        }
                    }
                    sumPerProc[i] = acc;
                    // timesPerProc.PushBack(timesOfThisProc, d.GetAllocator());
                    std::string indexStr = std::to_string(i);
                    rapidjson::Value indexStrRJ(indexStr.c_str(), indexStr.size(), d.GetAllocator());
                    timesPerProc.AddMember(indexStrRJ, timesOfThisProc, d.GetAllocator());
                }

                int64_t sumOfAllProc = 0;
                for (int64_t& t : sumPerProc) {
                    // check for over & underflow. https://stackoverflow.com/a/1514309
                    if (t > 0 && sumOfAllProc > std::numeric_limits<int64_t>::max() - t) {
                        std::cout << "Overflow Warning for profiling" << std::endl;
                    }
                    if (t < 0 && sumOfAllProc < std::numeric_limits<int64_t>::min() - t) {
                        std::cout << "Underflow Warning for profiling" << std::endl;
                    }
                    sumOfAllProc += t;
                }

                double sumOfAvgUnionProcTimes = 0;
                int64_t sumOfMaxUnionProcTimes = 0;
                for (int i = 0; i < maxNumElements; i++) {
                    int64_t sumSubStep = 0;
                    int64_t maxVal = std::numeric_limits<int64_t>::min();
                    int64_t minVal = std::numeric_limits<int64_t>::max();
                    int divider = 0;
                    for (int j = 0; j < simulation->GetTopology()->GetWorldSize(); j++) {
                        if (i < numAllElements[j]) {
                            int64_t t = allElements[displacements[j] + i];
                            sumSubStep += t;
                            divider++;

                            if (t > maxVal) {
                                maxVal = t;
                            }

                            if (t < minVal) {
                                minVal = t;
                            }
                        }
                    }
                    sumOfAvgUnionProcTimes += (double)sumSubStep / (double)divider;
                    sumOfMaxUnionProcTimes += maxVal;
                }

                // calculate values for stvr
                if (k.compare("CalculateForces") == 0 || k.compare("calculateForces") == 0) {
                    for (size_t i = 0; i < sumPerProc.size(); i++) {
                        double val = (double)sumPerProc[i] / (double)numAllElements[i];
                        valuesForSTVR.push_back(val);
                    }
                }

                // TODO: this is maybe incorrect for the NATA algorithm
                double avgTimeAllProc =
                    (double)sumOfAllProc / ((double)simulation->GetTopology()->GetWorldSize() * (double)maxNumElements);
                double totalTimeAllProc = (double)sumOfAllProc;

                rapidjson::Value oInner(rapidjson::kObjectType);

                rapidjson::Value timeUnitRJ(rapidjson::kStringType);
                rapidjson::Value avgTimeAllProcRJ(rapidjson::kNumberType);
                rapidjson::Value totalTimeAllProcRJ(rapidjson::kNumberType);
                rapidjson::Value totalAvgTimeUnionProcRJ(rapidjson::kNumberType);
                rapidjson::Value totalMaxTimeUnionProcRJ(rapidjson::kNumberType);

                std::string tu = charToTimeUnit(v.first);
                timeUnitRJ.SetString(tu.c_str(), tu.size(), d.GetAllocator());
                avgTimeAllProcRJ.SetDouble(avgTimeAllProc);
                totalTimeAllProcRJ.SetDouble(totalTimeAllProc);
                totalAvgTimeUnionProcRJ.SetDouble(sumOfAvgUnionProcTimes);
                totalMaxTimeUnionProcRJ.SetInt64(sumOfMaxUnionProcTimes);

                oInner.AddMember("time unit", timeUnitRJ, d.GetAllocator());
                oInner.AddMember("avg time all proc", avgTimeAllProcRJ, d.GetAllocator());
                oInner.AddMember("total time all proc", totalTimeAllProcRJ, d.GetAllocator());
                oInner.AddMember("total avg time union proc", totalAvgTimeUnionProcRJ, d.GetAllocator());
                oInner.AddMember("total max time union proc", totalMaxTimeUnionProcRJ, d.GetAllocator());
                oInner.AddMember("times per proc", timesPerProc, d.GetAllocator());

                rapidjson::Value key(k.c_str(), k.size(), d.GetAllocator());
                o.AddMember(key, oInner, d.GetAllocator());
            }
        }
    }

    // hitrate calculation
    {
        std::vector<double> allHitrates;
        int numMyHitrates = hitrates.size();
        std::vector<int> numAllHitrates;
        numAllHitrates.resize(simulation->GetTopology()->GetWorldSize());
        // elements from each process are gathered in order of their rank
        MPI_Gather(&numMyHitrates, 1, MPI_INT, numAllHitrates.data(), 1, MPI_INT, 0,
                   simulation->GetTopology()->GetComm());

        int maxNumElements = *(std::max_element(numAllHitrates.begin(), numAllHitrates.end()));

        std::vector<int> displacements;
        int sumDispl = 0;
        for (int i = 0; i < simulation->GetTopology()->GetWorldSize(); i++) {
            displacements.push_back(sumDispl);
            sumDispl += maxNumElements;
        }

        allHitrates.resize(sumDispl);
        std::fill(allHitrates.begin(), allHitrates.end(), 0.0);

        MPI_Gatherv(hitrates.data(), hitrates.size(), MPI_DOUBLE, allHitrates.data(), numAllHitrates.data(),
                    displacements.data(), MPI_DOUBLE, 0, simulation->GetTopology()->GetComm());

        if (simulation->GetTopology()->GetWorldRank() == 0) {
            /*for (size_t i = 0; i < numAllHitrates.size(); i++) {
                std::cout << "proc " << i << " has a numHitrates of " << numAllHitrates[i] << std::endl;
            }*/

            int sumNumAllHitrates = 0;
            for (int& nHr : numAllHitrates) {
                sumNumAllHitrates += nHr;
            }

            std::vector<double> accHitrates;
            accHitrates.resize(maxNumElements);

            rapidjson::Value hitratesPerSimStepRJ(rapidjson::kArrayType);

            for (int j = 0; j < maxNumElements; j++) {
                double accHrPerSimStep = 0;
                for (int i = 0; i < simulation->GetTopology()->GetWorldSize(); i++) {
                    // std::cout << "hr of proc " << i << " in step " << j << ": " << allHitrates[displacements[i] + j]
                    //          << std::endl;
                    double val = allHitrates[displacements[i] + j];
                    accHrPerSimStep += val;

                    if (j < numAllHitrates[i]) {
                        rapidjson::Value hrOfThisProcRJ(val);

                        std::string indexStr = std::to_string(i);
                        rapidjson::Value indexStrRJ(indexStr.c_str(), indexStr.size(), d.GetAllocator());
                        hrPerProcRJ.AddMember(indexStrRJ, hrOfThisProcRJ, d.GetAllocator());
                    }
                }
                // Note: this works only if we have just one sim step
                // TODO: make this work with multiple time steps
                accHitrates[j] = accHrPerSimStep / (double)sumNumAllHitrates;
                // std::cout << "divided hr by " << sumNumAllHitrates << std::endl;
            }

            // Note: this works only if we have just one sim step
            hitratesPerSimStepRJ.PushBack(rapidjson::Value(accHitrates[0]), d.GetAllocator());

            // rapidjson::Value hitrateRJ(rapidjson::kNumberType);
            // hitrateRJ.SetDouble(
            //    sum / ((double)simulation->GetTopology()->GetWorldSize() * (double)simulation->GetNumIterations()));
            // hr.AddMember("avg hitrate over all sim steps", hitrateRJ, d.GetAllocator());
            hr.AddMember("avg hitrate per sim step", hitratesPerSimStepRJ, d.GetAllocator());
        }
    }

    // STVR
    {
        if (simulation->GetTopology()->GetWorldRank() == 0) {
            double avgStepTimeAllProc = 0;
            for (const double& t_j : valuesForSTVR) {
                avgStepTimeAllProc += t_j;
            }
            avgStepTimeAllProc /= (double)(simulation->GetTopology()->GetWorldSize());

            int i = 0;
            for (const double& t_i : valuesForSTVR) {
                double stvr = std::abs((t_i - avgStepTimeAllProc) / avgStepTimeAllProc);

                // stvrRJ.PushBack(rapidjson::Value(stvr), d.GetAllocator());

                rapidjson::Value stvrPerProcRJ(rapidjson::kObjectType);
                rapidjson::Value cartRankRJ(rapidjson::kArrayType);

                CartRank rank = cartRanks[i];

                for (int j = 0; j < rank.GetDimensions(); j++) {
                    // std::cout << rank.GetRank()[j] << std::endl;
                    cartRankRJ.PushBack(rapidjson::Value(rank.GetRank()[j]), d.GetAllocator());
                }

                stvrPerProcRJ.AddMember("rank", cartRankRJ, d.GetAllocator());
                stvrPerProcRJ.AddMember("stvr", rapidjson::Value(stvr), d.GetAllocator());

                stvrRJ.PushBack(stvrPerProcRJ, d.GetAllocator());
                i++;
            }
        }
    }

    if (simulation->GetTopology()->GetWorldRank() == 0) {
        d.AddMember("times", o, d.GetAllocator());
        d.AddMember("hitrate", hr, d.GetAllocator());
        d.AddMember("hitrate per proc", hrPerProcRJ, d.GetAllocator());
        d.AddMember("stvr", stvrRJ, d.GetAllocator());

        rapidjson::StringBuffer buffer;
        rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);
        d.Accept(writer);

        std::ofstream csvFile;
        csvFile.open(outFile);
        csvFile << buffer.GetString();
        csvFile.close();
    }
}

#endif

#if defined(VLEVEL) && !defined(BENCHMARK_3BMDA) && !defined(TESTS_3BMDA)
void gatherAndPrintMessages()
{
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
#endif

int main(int argc, char* argv[])
{
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

    // decompositions = Utility::jsonToDecompositionArray("decompositions.json", "optimal");

    std::shared_ptr<Simulation> simulation;

    switch (a.algorithm) {
        case AlgorithmType::NATAType: simulation = createNATAContext(a.outputCSV); break;
        case AlgorithmType::P3BCAType:
            simulation = createP3BCAContext(
                a.outputCSV,
                Utility::getDecomposition(worldSize, (a.optimalDecomposition ? decompositions : decompositionsNaive)));
            break;
        case AlgorithmType::AUTAType: simulation = createAUTAContext(a.outputCSV); break;
        default: simulation = createNATAContext(a.outputCSV); break;
    }

    simulation->Init();

    // execute simulation
    simulation->Start();

#if defined(VLEVEL) && !defined(BENCHMARK_3BMDA) && !defined(TESTS_3BMDA) && VLEVEL > 0

    std::string message = "I'm proc " + std::to_string(simulation->GetTopology()->GetWorldRank()) + " and have done " +
                          std::to_string(simulation->GetNumParticleInteractions(0)) +
                          " particle interactions actually, and " +
                          std::to_string(simulation->GetNumBufferInteractions(0)) + " buffer interactions";

    MPIReporter::instance()->StoreMessage(simulation->GetTopology()->GetWorldRank(), message);
#endif

#ifdef PROFILE_3BMDA
    doTimingStuff(simulation, a.outputProfile);
#endif

#if defined(VLEVEL) && !defined(BENCHMARK_3BMDA) && !defined(TESTS_3BMDA)
    gatherAndPrintMessages();
#endif

    // finalize
    MPI_Type_free(&mpiParticleType);
    MPI_Finalize();

    return 0;
}