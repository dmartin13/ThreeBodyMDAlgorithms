#pragma once
#include <chrono>
#include <vector>

#include "Algorithm.hpp"
#include "topology/RingTopology.hpp"

class NATA final : public Algorithm {
private:
    int leftNeighbor;
    int rightNeighbor;
    int worldRank;
    int worldSize;
    int b1Owner;
    int b2Owner;

    std::shared_ptr<RingTopology> ringTopology;

    std::vector<Utility::Particle> b0;
    std::vector<Utility::Particle> b1;
    std::vector<Utility::Particle> b2;

    std::vector<Utility::Triplet> alreadyProcessed;

    bool containsProcessed(Utility::Triplet t);
    void calculateProcessed(int step, bool &calculate);
    int shiftRight(std::vector<Utility::Particle> &buf);

public:
    NATA();
    virtual ~NATA();

    void Init(std::shared_ptr<Simulation> simulation) override;

    std::tuple<uint64_t, uint64_t> SimulationStep() override;
    std::tuple<uint64_t, uint64_t> SimulationStep([[maybe_unused]] ForceType forceType) override;
};
