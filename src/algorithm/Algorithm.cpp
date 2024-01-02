#include "Algorithm.hpp"

Algorithm::Algorithm() {}

Algorithm::~Algorithm() {}

void Algorithm::Init(std::shared_ptr<Simulation> simulation) {
    this->simulation = simulation;
    this->mpiParticleType = simulation->GetMPIParticleType();
    this->potential = this->simulation->GetTriwisePotential();
    this->pairwisePotential = this->simulation->GetPairwisePotential();
    this->worldSize = this->simulation->GetTopology()->GetWorldSize();
    this->worldRank = this->simulation->GetTopology()->GetWorldRank();
}

std::tuple<uint64_t, uint64_t> Algorithm::CalculateInteractions(std::vector<Utility::Particle> &b0,
                                                                std::vector<Utility::Particle> &b1,
                                                                std::vector<Utility::Particle> &b2, int b0Owner,
                                                                int b1Owner, int b2Owner) {
    return calculateInteractions(b0, b1, b2, b0Owner, b1Owner, b2Owner, 0, -1);
}

std::tuple<uint64_t, uint64_t> Algorithm::CalculateInteractions(std::vector<Utility::Particle> &b0,
                                                                std::vector<Utility::Particle> &b1,
                                                                std::vector<Utility::Particle> &b2, int b0Owner,
                                                                int b1Owner, int b2Owner, int b0Start, int b0NumSteps) {
    return calculateInteractions(b0, b1, b2, b0Owner, b1Owner, b2Owner, b0Start, b0NumSteps);
}

std::tuple<uint64_t, uint64_t> Algorithm::calculateInteractions(std::vector<Utility::Particle> &b0,
                                                                std::vector<Utility::Particle> &b1,
                                                                std::vector<Utility::Particle> &b2, int b0Owner,
                                                                int b1Owner, int b2Owner, int b0Start, int b0NumSteps) {
    uint64_t numActParticleInteractions = 0;
    uint64_t numPossibleParticleInteractions = 0;
    for (size_t i = b0Start; i < (b0NumSteps != -1 ? (size_t)(b0Start + b0NumSteps) : b0.size()); ++i) {
        if (b0[i].isDummy) {
            continue;
        }
        int b1LoopIndex = b1Owner == b0Owner ? i + 1 : 0;
        for (size_t j = b1LoopIndex; j < b1.size(); ++j) {
            if (b1[j].isDummy) {
                continue;
            }
            int b2LoopIndex = 0;
            if (b2Owner == b1Owner) {
                b2LoopIndex = j + 1;
            } else if (b2Owner == b0Owner) {
                b2LoopIndex = i + 1;
            }
            for (size_t k = b2LoopIndex; k < b2.size(); ++k) {
                if (b2[k].isDummy) {
                    continue;
                }
                numPossibleParticleInteractions++;

                potential->CalculateForces(b0[i], b1[j], b2[k]);

                numActParticleInteractions++;
            }
        }
    }

    return std::tuple(numActParticleInteractions, numPossibleParticleInteractions);
}

void Algorithm::SumUpParticles(std::vector<Utility::Particle> &b0, std::vector<Utility::Particle> &b1,
                               std::vector<Utility::Particle> &b2) {
    for (size_t i = 0; i < b0.size(); i++) {
        b0[i].f0X += b1[i].f0X + b2[i].f0X;
        b0[i].f0Y += b1[i].f0Y + b2[i].f0Y;
        b0[i].f0Z += b1[i].f0Z + b2[i].f0Z;

        b0[i].f1X += b1[i].f1X + b2[i].f1X;
        b0[i].f1Y += b1[i].f1Y + b2[i].f1Y;
        b0[i].f1Z += b1[i].f1Z + b2[i].f1Z;
    }
}

void Algorithm::SumUpParticles(std::vector<Utility::Particle> &b0, std::vector<Utility::Particle> &b1) {
    for (size_t i = 0; i < b0.size(); i++) {
        b0[i].f0X += b1[i].f0X;
        b0[i].f0Y += b1[i].f0Y;
        b0[i].f0Z += b1[i].f0Z;

        b0[i].f1X += b1[i].f1X;
        b0[i].f1Y += b1[i].f1Y;
        b0[i].f1Z += b1[i].f1Z;
    }
}

std::tuple<uint64_t, uint64_t> Algorithm::CalculatePairwiseInteractions(std::vector<Utility::Particle> &b0,
                                                                        std::vector<Utility::Particle> &b1, int b0Owner,
                                                                        int b1Owner) {
    return calculatePairwiseInteractions(b0, b1, b0Owner, b1Owner, 0, -1);
}

std::tuple<uint64_t, uint64_t> Algorithm::calculatePairwiseInteractions(std::vector<Utility::Particle> &b0,
                                                                        std::vector<Utility::Particle> &b1, int b0Owner,
                                                                        int b1Owner, int b0Start, int b0NumSteps) {
    uint64_t numActParticleInteractions = 0;
    uint64_t numPossibleParticleInteractions = 0;
    for (size_t i = b0Start; i < (b0NumSteps != -1 ? (size_t)(b0Start + b0NumSteps) : b0.size()); ++i) {
        if (b0[i].isDummy) {
            continue;
        }
        int b1LoopIndex = b1Owner == b0Owner ? i + 1 : 0;
        for (size_t j = b1LoopIndex; j < b1.size(); ++j) {
            if (b1[j].isDummy) {
                continue;
            }
            pairwisePotential->CalculateForces(b0[i], b1[j]);
        }
    }

    return std::tuple(numActParticleInteractions, numPossibleParticleInteractions);
}

int Algorithm::GetNumShifts() { return this->numShifts; }