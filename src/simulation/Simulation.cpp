#include "Simulation.hpp"

Simulation::Simulation(int iterations, std::shared_ptr<Algorithm> algorithm, std::shared_ptr<Topology> topology,
                       std::shared_ptr<PairwisePotential> pairwisepotential,
                       std::shared_ptr<TriwisePotential> triwisepotential,
                       std::shared_ptr<DomainDecomposition> decomposition, MPI_Datatype* mpiParticleType,
                       std::vector<Utility::Particle>& particles, double dt, Eigen::Vector3d gForce,
                       std::string csvOutput)
    : iterations(iterations), algorithm(algorithm), topology(topology), pairwisepotential(pairwisepotential),
      triwisepotential(triwisepotential), decomposition(decomposition), mpiParticleType(mpiParticleType),
      particles(particles), dt(dt), gForce(gForce), csvOutput(csvOutput)
{}

void Simulation::Init()
{
    std::shared_ptr<Simulation> simulationPtr = shared_from_this();
    this->topology->Init(simulationPtr);
    this->decomposition->Init(simulationPtr);
    this->algorithm->Init(simulationPtr);
    // this->potential->Init(simulationPtr);
    this->pairwisepotential->Init(simulationPtr);
    this->triwisepotential->Init(simulationPtr);
}

Simulation::~Simulation() {}

std::shared_ptr<Algorithm> Simulation::GetAlgorithm() { return this->algorithm; }
std::shared_ptr<Topology> Simulation::GetTopology() { return this->topology; }
std::shared_ptr<PairwisePotential> Simulation::GetPairwisePotential() { return this->pairwisepotential; }
std::shared_ptr<TriwisePotential> Simulation::GetTriwisePotential() { return this->triwisepotential; }
std::shared_ptr<DomainDecomposition> Simulation::GetDecomposition() { return this->decomposition; }

void Simulation::Start()
{
    for (int i = 0; i < iterations; ++i) {
#ifdef MEASURESIMSTEP_3BMDA
        std::chrono::time_point<std::chrono::system_clock> start;
        std::chrono::time_point<std::chrono::system_clock> end;
        start = std::chrono::system_clock::now();
#endif

        // update particle positions at predictor stage
        this->decomposition->UpdatePredictorStage(this->dt);
        MPI_Barrier(this->topology->GetComm());

        // execute algorithm... force calculation
        numInteractions.push_back(this->algorithm->SimulationStep());
        MPI_Barrier(this->topology->GetComm());

        // update the particle positions
        this->decomposition->Update(this->dt, this->gForce);
        MPI_Barrier(this->topology->GetComm());

#ifdef MEASURESIMSTEP_3BMDA
        end = std::chrono::system_clock::now();
        auto elapsed_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
        std::vector<uint64_t> allTimes;
        if (topology->GetWorldRank() == 0) {
            allTimes.resize(topology->GetWorldSize());
        }
        MPI_Gather(&elapsed_time, 1, MPI_INT64_T, allTimes.data(), 1, MPI_INT64_T, 0, topology->GetComm());
        if (topology->GetWorldRank() == 0) {
            uint64_t max = *max_element(allTimes.begin(), allTimes.end());
            std::string allTimesStr = "[";
            for (size_t i = 0; i < allTimes.size(); i++) {
                allTimesStr.append(std::to_string(allTimes[i]));
                if (i < allTimes.size() - 1) {
                    allTimesStr.append(", ");
                }
            }
            allTimesStr.append("]");
            std::string json = "{\"allTimes\": " + allTimesStr + ", \"max\": " + std::to_string(max) + "}";
            std::cout << json << std::endl;
        }

#endif

#if !defined(BENCHMARK_3BMDA) && !defined(TESTS_3BMDA) && !defined(PROFILE_3BMDA)
        if (this->csvOutput.compare("") != 0) {
            writeSimulationStepToCSV(this->csvOutput.substr(0, this->csvOutput.find_last_of('.')) + "_" +
                                     std::to_string(i) + ".csv");
        }
#endif
    }
}

MPI_Datatype* Simulation::GetMPIParticleType() { return this->mpiParticleType; }

std::vector<Utility::Particle>& Simulation::GetAllParticles() { return this->particles; }

double Simulation::GetDeltaT() { return this->dt; }
int Simulation::GetNumIterations() { return this->iterations; }
uint64_t Simulation::GetNumBufferInteractions(int step)
{
    return (size_t)step < this->numInteractions.size() ? std::get<0>(this->numInteractions[step]) : 0;
}
uint64_t Simulation::GetNumParticleInteractions(int step)
{
    return (size_t)step < this->numInteractions.size() ? std::get<1>(this->numInteractions[step]) : 0;
}
Eigen::Vector3d Simulation::GetGForce() { return this->gForce; }

void Simulation::writeSimulationStepToCSV(std::string file)
{
    std::vector<Utility::Particle> receivedParticlesForCSVOutput;

    int numProcessors = topology->GetWorldSize();
    std::vector<int> numParticlesPerProcessor;
    numParticlesPerProcessor.resize(numProcessors);
    // elements from each process are gathered in order of their rank
    int numOfMyParticles = decomposition->GetNumOfMyParticles();
    MPI_Gather(&numOfMyParticles, 1, MPI_INT, numParticlesPerProcessor.data(), 1, MPI_INT, 0, topology->GetComm());

    std::vector<int> displacements;
    int sumDispl = 0;
    for (int i = 0; i < numProcessors; i++) {
        displacements.push_back(sumDispl);
        sumDispl += numParticlesPerProcessor[i];
    }

    if (topology->GetWorldRank() == 0) {
        receivedParticlesForCSVOutput.resize(sumDispl);
    }

    std::vector<Utility::Particle> particlesToSend = this->decomposition->GetMyParticles();

    MPI_Gatherv(particlesToSend.data(), particlesToSend.size(), *mpiParticleType, receivedParticlesForCSVOutput.data(),
                numParticlesPerProcessor.data(), displacements.data(), *mpiParticleType, 0, topology->GetComm());

    std::sort(receivedParticlesForCSVOutput.begin(), receivedParticlesForCSVOutput.end(),
              [](Utility::Particle a, Utility::Particle b) { return a.ID < b.ID; });

    if (topology->GetWorldRank() == 0) {
        Utility::writeStepToCSVWithForces(file, receivedParticlesForCSVOutput);
    }
}