#pragma once

#ifdef MEASURESIMSTEP_3BMDA
#include <chrono>
#endif

#include "algorithm/Algorithm.hpp"
#include "fwd.hpp"
#include "potential/PairwisePotential.hpp"
#include "potential/Potential.hpp"
#include "potential/TriwisePotential.hpp"

class Simulation : public std::enable_shared_from_this<Simulation> {
private:
    int iterations;
    std::shared_ptr<Algorithm> algorithm;
    std::shared_ptr<Topology> topology;
    // std::shared_ptr<Potential> potential;
    std::shared_ptr<PairwisePotential> pairwisepotential;
    std::shared_ptr<TriwisePotential> triwisepotential;
    std::shared_ptr<DomainDecomposition> decomposition;
    MPI_Datatype* mpiParticleType;
    std::vector<Utility::Particle>& particles;
    double dt;
    Eigen::Vector3d gForce;
    std::vector<std::tuple<uint64_t, uint64_t>> numInteractions;
    std::string csvOutput;

    void writeSimulationStepToCSV(std::string file);

public:
    Simulation(int iterations, std::shared_ptr<Algorithm> algorithm, std::shared_ptr<Topology> topology,
               std::shared_ptr<PairwisePotential> pairwisepotential, std::shared_ptr<TriwisePotential> triwisepotential,
               std::shared_ptr<DomainDecomposition> decomposition, MPI_Datatype* mpiParticleType,
               std::vector<Utility::Particle>& particles, double dt, Eigen::Vector3d gForce,
               std::string csvOutput = "");
    virtual ~Simulation();

    void Start();
    void Init();

    std::shared_ptr<Algorithm> GetAlgorithm();
    std::shared_ptr<Topology> GetTopology();
    // std::shared_ptr<Potential> GetPotential();
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