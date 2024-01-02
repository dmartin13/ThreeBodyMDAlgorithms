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

    void calculateTriwiseInteractions(std::vector<Utility::Particle> &b0, std::vector<Utility::Particle> &b1,
                               std::vector<Utility::Particle> &b2, int b0Owner, int b1Owner, int b2Owner, int b0Start,
                               int b0NumSteps);

    void calculatePairwiseInteractions(std::vector<Utility::Particle> &b0, std::vector<Utility::Particle> &b1,
                                       int b0Owner, int b1Owner, int b0Start, int b0NumSteps);

public:
    Algorithm();
    virtual ~Algorithm();

    virtual void Init(std::shared_ptr<Simulation> simulation);

    virtual void SimulationStep() = 0;

    virtual void SimulationStep(ForceType forceType) = 0;

    void CalculateInteractions(std::vector<Utility::Particle> &b0, std::vector<Utility::Particle> &b1,
                               std::vector<Utility::Particle> &b2, int b0Owner, int b1Owner, int b2Owner);
    void CalculateInteractions(std::vector<Utility::Particle> &b0, std::vector<Utility::Particle> &b1,
                               std::vector<Utility::Particle> &b2, int b0Owner, int b1Owner, int b2Owner, int b0Start,
                               int b0NumSteps);

    void CalculatePairwiseInteractions(std::vector<Utility::Particle> &b0, std::vector<Utility::Particle> &b1,
                                       int b0Owner, int b1Owner);

    void SumUpParticles(std::vector<Utility::Particle> &b0, std::vector<Utility::Particle> &b1,
                        std::vector<Utility::Particle> &b2);

    void SumUpParticles(std::vector<Utility::Particle> &b0, std::vector<Utility::Particle> &b1);
};
