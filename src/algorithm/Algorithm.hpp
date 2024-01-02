#pragma once

#include <memory>
#include <utility/enums.hpp>

#include "MPIReporter.hpp"
#include "decomposition/DomainDecomposition.hpp"
#include "fwd.hpp"
#include "potential/PairwisePotential.hpp"
#include "potential/Potential.hpp"
#include "potential/TriwisePotential.hpp"
#include "simulation/Simulation.hpp"
#include "topology/Topology.hpp"
#include "utility/utility.hpp"

class Algorithm {
protected:
    std::shared_ptr<Simulation> simulation;
    MPI_Datatype *mpiParticleType;
    std::shared_ptr<TriwisePotential> potential;
    std::shared_ptr<PairwisePotential> pairwisePotential;
    int worldSize;
    int worldRank;
    int numShifts;

    std::tuple<uint64_t, uint64_t> calculateInteractions(std::vector<Utility::Particle> &b0,
                                                         std::vector<Utility::Particle> &b1,
                                                         std::vector<Utility::Particle> &b2, int b0Owner, int b1Owner,
                                                         int b2Owner, int b0Start, int b0NumSteps);

    std::tuple<uint64_t, uint64_t> calculatePairwiseInteractions(std::vector<Utility::Particle> &b0,
                                                                 std::vector<Utility::Particle> &b1, int b0Owner,
                                                                 int b1Owner, int b0Start, int b0NumSteps);

    void calcParticleInteractions(std::vector<std::tuple<int, int, int>> &particleTripletsToCalculate,
                                  std::vector<Utility::Particle> &b0, std::vector<Utility::Particle> &b1,
                                  std::vector<Utility::Particle> &b2);

public:
    Algorithm();
    virtual ~Algorithm();

    virtual void Init(std::shared_ptr<Simulation> simulation);

    virtual std::tuple<uint64_t, uint64_t> SimulationStep() = 0;

    virtual std::tuple<uint64_t, uint64_t> SimulationStep(ForceType forceType) = 0;

    std::tuple<uint64_t, uint64_t> CalculateInteractions(std::vector<Utility::Particle> &b0,
                                                         std::vector<Utility::Particle> &b1,
                                                         std::vector<Utility::Particle> &b2, int b0Owner, int b1Owner,
                                                         int b2Owner);
    std::tuple<uint64_t, uint64_t> CalculateInteractions(std::vector<Utility::Particle> &b0,
                                                         std::vector<Utility::Particle> &b1,
                                                         std::vector<Utility::Particle> &b2, int b0Owner, int b1Owner,
                                                         int b2Owner, int b0Start, int b0NumSteps);

    std::tuple<uint64_t, uint64_t> CalculatePairwiseInteractions(std::vector<Utility::Particle> &b0,
                                                                 std::vector<Utility::Particle> &b1, int b0Owner,
                                                                 int b1Owner);

    void SumUpParticles(std::vector<Utility::Particle> &b0, std::vector<Utility::Particle> &b1,
                        std::vector<Utility::Particle> &b2);

    void SumUpParticles(std::vector<Utility::Particle> &b0, std::vector<Utility::Particle> &b1);

    int GetNumShifts();
};
