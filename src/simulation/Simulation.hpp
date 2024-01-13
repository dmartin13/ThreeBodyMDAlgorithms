#pragma once

#ifdef MEASURESIMSTEP_3BMDA
#include <chrono>
#endif

#include <format>

#include "Thermostat.hpp"
#include "algorithm/Algorithm.hpp"
#include "fwd.hpp"
#include "potential/PairwisePotential.hpp"
#include "potential/Potential.hpp"
#include "potential/TriwisePotential.hpp"
#include "utility/cli.hpp"

class Simulation : public std::enable_shared_from_this<Simulation> {
private:
    Utility::cliArguments args;
    int iterations;
    size_t respaStepSize;
    std::shared_ptr<Algorithm> algorithm;
    std::shared_ptr<Topology> topology;
    std::shared_ptr<PairwisePotential> pairwisepotential;
    std::shared_ptr<TriwisePotential> triwisepotential;
    std::shared_ptr<DomainDecomposition> decomposition;
    MPI_Datatype* mpiParticleType;
    std::vector<Utility::Particle>& particles;
    double dt;
    Eigen::Vector3d gForce;
    std::vector<std::tuple<uint64_t, uint64_t>> numInteractions;
    std::string csvOutput;

    std::vector<double> kineticEnergy;
    std::vector<double> potentialEnergy;
    std::vector<double> totalEnergy;

    void writeSimulationStepToCSV(std::string file);
    double calculateKineticEnergy();

public:
    Simulation(Utility::cliArguments& args, int iterations, int respaStepSize, std::shared_ptr<Algorithm> algorithm,
               std::shared_ptr<Topology> topology, std::shared_ptr<PairwisePotential> pairwisepotential,
               std::shared_ptr<TriwisePotential> triwisepotential, std::shared_ptr<DomainDecomposition> decomposition,
               MPI_Datatype* mpiParticleType, std::vector<Utility::Particle>& particles, double dt,
               Eigen::Vector3d gForce, std::string csvOutput = "");
    virtual ~Simulation();

    void Start();
    void Init();

    std::shared_ptr<Algorithm> GetAlgorithm();
    std::shared_ptr<Topology> GetTopology();
    std::shared_ptr<PairwisePotential> GetPairwisePotential();
    std::shared_ptr<TriwisePotential> GetTriwisePotential();
    std::shared_ptr<DomainDecomposition> GetDecomposition();
    MPI_Datatype* GetMPIParticleType();
    std::vector<Utility::Particle>& GetAllParticles();
    double GetDeltaT();
    int GetNumIterations();
    uint64_t GetNumBufferInteractions(int step);
    uint64_t GetNumParticleInteractions(int step);

    Eigen::Vector3d GetGForce();
};