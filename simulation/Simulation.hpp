#pragma once

#include "../algorithm/Algorithm.hpp"
#include "../fwd.hpp"
#include "../potential/Potential.hpp"

class Simulation : public std::enable_shared_from_this<Simulation> {
private:
    int iterations;
    std::shared_ptr<Algorithm> algorithm;
    std::shared_ptr<Topology> topology;
    std::shared_ptr<Potential> potential;
    std::shared_ptr<DomainDecomposition> decomposition;
    MPI_Datatype* mpiParticleType;
    std::vector<Utility::Particle>& particles;
    double dt;
    Eigen::Vector3d gForce;
    std::vector<std::tuple<int, int>> numInteractions;
    std::string csvOutput;

    void writeSimulationStepToCSV(std::string file);

#ifdef PROFILE_3BMDA
    //void calcOverallHitrate(int step, double hitrate);
#endif

public:
    Simulation(int iterations, std::shared_ptr<Algorithm> algorithm, std::shared_ptr<Topology> topology,
               std::shared_ptr<Potential> potential, std::shared_ptr<DomainDecomposition> decomposition,
               MPI_Datatype* mpiParticleType, std::vector<Utility::Particle>& particles, double dt,
               Eigen::Vector3d gForce, std::string csvOutput = "");
    virtual ~Simulation();

    void Start();
    void Init();

    std::shared_ptr<Algorithm> GetAlgorithm();
    std::shared_ptr<Topology> GetTopology();
    std::shared_ptr<Potential> GetPotential();
    std::shared_ptr<DomainDecomposition> GetDecomposition();
    MPI_Datatype* GetMPIParticleType();
    std::vector<Utility::Particle>& GetAllParticles();
    double GetDeltaT();
    int GetNumIterations();
    int GetNumBufferInteractions(int step);
    int GetNumParticleInteractions(int step);

    Eigen::Vector3d GetGForce();
};