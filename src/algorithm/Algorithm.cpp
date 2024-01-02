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
#ifdef PROFILE_3BMDA
    this->calcForcesAcc = 0;
#endif

#if defined(VLEVEL) && !defined(BENCHMARK_3BMDA) && !defined(TESTS_3BMDA) && VLEVEL > 0
    std::string message = "I'm proc " + std::to_string(this->simulation->GetTopology()->GetWorldRank()) + " and own " +
                          std::to_string(this->simulation->GetDecomposition()->GetMyParticles().size()) + " particles";
    MPIReporter::instance()->StoreMessage(this->simulation->GetTopology()->GetWorldRank(), message);

#endif
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
#ifdef PROFILE_3BMDA
    this->calcForcesAcc = 0;
    // bool append = false;
    std::chrono::time_point<std::chrono::system_clock> start;
    std::chrono::time_point<std::chrono::system_clock> end;
    start = std::chrono::system_clock::now();
#endif
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

#ifdef PROFILE_3BMDA
                std::chrono::time_point<std::chrono::system_clock> start1;
                std::chrono::time_point<std::chrono::system_clock> end1;
                start1 = std::chrono::system_clock::now();
#endif
                potential->CalculateForces(b0[i], b1[j], b2[k]);
#ifdef PROFILE_3BMDA
                end1 = std::chrono::system_clock::now();
                auto elapsed_time1 = std::chrono::duration_cast<std::chrono::nanoseconds>(end1 - start1);
                this->calcForcesAcc += elapsed_time1.count();
#endif
                numActParticleInteractions++;
            }
        }
    }
#ifdef PROFILE_3BMDA
    end = std::chrono::system_clock::now();
    auto elapsed_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    // std::cout << "proc " << this->simulation->GetTopology()->GetWorldRank() << " calcstep took " <<
    // elapsed_time.count()
    //          << " ns" << std::endl;
    bool hasKey = this->times.count("calculateInteractions");
    if (!hasKey) {
        this->times["calculateInteractions"] = std::make_pair(0, std::vector<int64_t>());
    }
    this->times["calculateInteractions"].second.push_back(elapsed_time.count());

    bool hasKey2 = this->times.count("calculateForces");
    if (!hasKey2) {
        this->times["calculateForces"] = std::make_pair(0, std::vector<int64_t>());
    }
    this->times["calculateForces"].second.push_back(this->calcForcesAcc);
#endif

    return std::tuple(numActParticleInteractions, numPossibleParticleInteractions);
}

#ifdef PROFILE_3BMDA
void Algorithm::calcParticleInteractions(std::vector<std::tuple<int, int, int>> &particleTripletsToCalculate,
                                         std::vector<Utility::Particle> &b0, std::vector<Utility::Particle> &b1,
                                         std::vector<Utility::Particle> &b2, bool append)
#else
void Algorithm::calcParticleInteractions(std::vector<std::tuple<int, int, int>> &particleTripletsToCalculate,
                                         std::vector<Utility::Particle> &b0, std::vector<Utility::Particle> &b1,
                                         std::vector<Utility::Particle> &b2)
#endif
{
#ifdef PROFILE_3BMDA
    std::chrono::time_point<std::chrono::system_clock> start1;
    std::chrono::time_point<std::chrono::system_clock> end1;
    start1 = std::chrono::system_clock::now();
#endif

#if defined(USE_OMP) && defined(OPENMPAVAIL)
#pragma omp parallel for
    for (auto it = particleTripletsToCalculate.begin(); it < particleTripletsToCalculate.end(); it++) {
        this->potential->CalculateForces(b0[std::get<0>((*it))], b1[std::get<1>((*it))], b2[std::get<2>((*it))]);
    }
#else
    for (auto it = particleTripletsToCalculate.begin(); it < particleTripletsToCalculate.end(); it++) {
        potential->CalculateForces(b0[std::get<0>((*it))], b1[std::get<1>((*it))], b2[std::get<2>((*it))]);
    }
#endif

#ifdef PROFILE_3BMDA
    end1 = std::chrono::system_clock::now();
    auto elapsed_time1 = std::chrono::duration_cast<std::chrono::nanoseconds>(end1 - start1);
    bool hasKey = this->times.count("CalculateForces");
    if (!hasKey) {
        this->times["CalculateForces"] = std::make_pair(0, std::vector<int64_t>());
    }
    if (append && this->times["CalculateForces"].second.size() > 0) {
        this->times["CalculateForces"].second.back() += elapsed_time1.count();
    } else {
        this->times["CalculateForces"].second.push_back(elapsed_time1.count());
    }
#endif
}

void Algorithm::SumUpParticles(std::vector<Utility::Particle> &b0, std::vector<Utility::Particle> &b1,
                               std::vector<Utility::Particle> &b2) {
#ifdef PROFILE_3BMDA
    std::chrono::time_point<std::chrono::system_clock> start;
    std::chrono::time_point<std::chrono::system_clock> end;
    start = std::chrono::system_clock::now();
#endif
    for (size_t i = 0; i < b0.size(); i++) {
        b0[i].f0X += b1[i].f0X + b2[i].f0X;
        b0[i].f0Y += b1[i].f0Y + b2[i].f0Y;
        b0[i].f0Z += b1[i].f0Z + b2[i].f0Z;

        b0[i].f1X += b1[i].f1X + b2[i].f1X;
        b0[i].f1Y += b1[i].f1Y + b2[i].f1Y;
        b0[i].f1Z += b1[i].f1Z + b2[i].f1Z;
    }
#ifdef PROFILE_3BMDA
    end = std::chrono::system_clock::now();
    auto elapsed_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    bool hasKey = this->times.count("SumUpParticles");
    if (!hasKey) {
        this->times["SumUpParticles"] = std::make_pair(0, std::vector<int64_t>());
    }
    this->times["SumUpParticles"].second.push_back(elapsed_time.count());
#endif
}

void Algorithm::SumUpParticles(std::vector<Utility::Particle> &b0, std::vector<Utility::Particle> &b1) {
#ifdef PROFILE_3BMDA
    std::chrono::time_point<std::chrono::system_clock> start;
    std::chrono::time_point<std::chrono::system_clock> end;
    start = std::chrono::system_clock::now();
#endif
    for (size_t i = 0; i < b0.size(); i++) {
        b0[i].f0X += b1[i].f0X;
        b0[i].f0Y += b1[i].f0Y;
        b0[i].f0Z += b1[i].f0Z;

        b0[i].f1X += b1[i].f1X;
        b0[i].f1Y += b1[i].f1Y;
        b0[i].f1Z += b1[i].f1Z;
    }
#ifdef PROFILE_3BMDA
    end = std::chrono::system_clock::now();
    auto elapsed_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    bool hasKey = this->times.count("SumUpParticles");
    if (!hasKey) {
        this->times["SumUpParticles"] = std::make_pair(0, std::vector<int64_t>());
    }
    this->times["SumUpParticles"].second.push_back(elapsed_time.count());
#endif
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

#ifdef TESTS_3BMDA
std::vector<Utility::Triplet> Algorithm::GetProcessed() { return this->processed; }
#endif

#ifdef PROFILE_3BMDA
std::map<std::string, std::pair<char, std::vector<int64_t>>> Algorithm::GetTimes() { return this->times; }

std::vector<double> Algorithm::GetHitrates() { return this->hitrates; }
#endif

int Algorithm::GetNumShifts() { return this->numShifts; }