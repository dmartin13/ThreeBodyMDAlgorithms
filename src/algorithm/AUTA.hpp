#pragma once
#include <vector>

#include "topology/RingTopology.hpp"
#include "Algorithm.hpp"

class AUTA final : public Algorithm {
private:
    int leftNeighbor;
    int rightNeighbor;
    int worldRank;
    int worldSize;

    std::shared_ptr<RingTopology> ringTopology;

    std::vector<Utility::Particle> b0;
    std::vector<Utility::Particle> b1;
    std::vector<Utility::Particle> b2;
    std::vector<Utility::Particle> b0Tmp;
    std::vector<Utility::Particle> b1Tmp;
    std::vector<Utility::Particle> b2Tmp;
    int b0Owner;
    int b1Owner;
    int b2Owner;

    int shiftRight(std::vector<Utility::Particle>& buf, int owner);
    std::tuple<uint64_t, uint64_t> calculateOneThirdOfInteractions(int thirdID);
    std::vector<Utility::Particle>* pickBuffer(int i);
    int& getBufOwner(int i);
    void sendBackParticles();

public:
    AUTA();
    virtual ~AUTA();

    void Init(std::shared_ptr<Simulation> simulation) override;

    std::tuple<uint64_t, uint64_t> SimulationStep() override;
};
