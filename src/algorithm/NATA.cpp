#include "NATA.hpp"

NATA::NATA() {}

NATA::~NATA() {}

void NATA::Init(std::shared_ptr<Simulation> simulation) {
    Algorithm::Init(simulation);

    this->b0 = this->simulation->GetDecomposition()->GetMyParticles();
    this->ringTopology = (std::static_pointer_cast<RingTopology>(this->simulation->GetTopology()));
    this->leftNeighbor = ringTopology->GetLeftNeighbor();
    this->rightNeighbor = ringTopology->GetRightNeighbor();
    this->worldRank = ringTopology->GetWorldRank();
    this->worldSize = ringTopology->GetWorldSize();
}

bool NATA::containsProcessed(Utility::Triplet t) {
    if (std::find(alreadyProcessed.begin(), alreadyProcessed.end(), t) != alreadyProcessed.end()) {
        return true;
    }
    return false;
}

void NATA::calculateProcessed(int step, bool &calculate) {
#ifdef PROFILE_3BMDA
    std::chrono::time_point<std::chrono::system_clock> start;
    std::chrono::time_point<std::chrono::system_clock> end;
    start = std::chrono::system_clock::now();
#endif
    for (int i = 0; i < this->worldSize; i++) {
        int b1Rank = Utility::mod(i - (step / this->worldSize), this->worldSize);
        int b2Rank = Utility::mod(i - step, this->worldSize);
        Utility::Triplet t(i, b1Rank, b2Rank);

        if (!containsProcessed(t)) {
            alreadyProcessed.push_back(t);
            if (i == this->worldRank) {
                calculate = true;
            }
        }
    }
#ifdef PROFILE_3BMDA
    end = std::chrono::system_clock::now();
    auto elapsed_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    bool hasKey = this->times.count("calculateProcessed");
    if (!hasKey) {
        this->times["calculateProcessed"] = std::make_pair(0, std::vector<int64_t>());
    }
    this->times["calculateProcessed"].second.push_back(elapsed_time.count());
#endif
}

int NATA::shiftRight(std::vector<Utility::Particle> &buf) {
#ifdef PROFILE_3BMDA
    std::chrono::time_point<std::chrono::system_clock> start;
    std::chrono::time_point<std::chrono::system_clock> end;
    start = std::chrono::system_clock::now();
#endif
    // https://moodle.rrze.uni-erlangen.de/pluginfile.php/13157/mod_resource/content/1/06_MPI_Advanced.pdf Page 12
    MPI_Sendrecv_replace(buf.data(), buf.size(), *this->mpiParticleType, this->rightNeighbor, 0, this->leftNeighbor, 0,
                         this->ringTopology->GetComm(), MPI_STATUS_IGNORE);
#ifdef PROFILE_3BMDA
    end = std::chrono::system_clock::now();
    auto elapsed_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    bool hasKey = this->times.count("shiftRight");
    if (!hasKey) {
        this->times["shiftRight"] = std::make_pair(0, std::vector<int64_t>());
    }
    this->times["shiftRight"].second.push_back(elapsed_time.count());
#endif

    this->numShifts++;

    return 0;
}

std::tuple<uint64_t, uint64_t> NATA::SimulationStep() {
    this->b0 = this->simulation->GetDecomposition()->GetMyParticles();

    this->numShifts = 0;

    // use assignment operator to copy vector
    b1 = b0;
    b2 = b0;

    this->b1Owner = this->worldRank;
    this->b2Owner = this->worldRank;

#ifdef TESTS_3BMDA
    processed.clear();
#endif
#ifdef PROFILE_3BMDA
    this->hitrate = 0;
    uint64_t hitRateDivider = 0;
#endif

    uint64_t numBufferInteractions = 0;
    uint64_t numParticleInteractionsAcc = 0;

    bool calculate = false;
    alreadyProcessed.clear();
    int step = 0;
    for (int i = 0; i < this->worldSize; i++) {
        this->b1Owner = Utility::mod(this->worldRank - (step / this->worldSize), this->worldSize);
        for (int j = 0; j < this->worldSize; j++) {
            this->b2Owner = Utility::mod(this->worldRank - step, this->worldSize);

            calculate = false;
            calculateProcessed(step, calculate);
            if (calculate) {
#if defined(VLEVEL) && !defined(BENCHMARK_3BMDA) && !defined(TESTS_3BMDA) && VLEVEL > 0
                std::string message = "I'm proc " + std::to_string(simulation->GetTopology()->GetWorldRank()) +
                                      " and going to calculate interactions between buffers from proc (" +
                                      std::to_string(worldRank) + ", " + std::to_string(this->b1Owner) + ", " +
                                      std::to_string(this->b2Owner) + ")";
                MPIReporter::instance()->StoreMessage(this->simulation->GetTopology()->GetWorldRank(), message);
#endif

                std::tuple<uint64_t, uint64_t> numParticleInteractions = this->CalculateInteractions(
                    this->b0, this->b1, this->b2, this->worldRank, this->b1Owner, this->b2Owner);

                numParticleInteractionsAcc += std::get<0>(numParticleInteractions);
                numBufferInteractions++;
#ifdef TESTS_3BMDA
                processed.push_back(Utility::Triplet(this->worldRank,
                                                     Utility::mod(i + this->worldRank, this->worldSize),
                                                     Utility::mod(j + this->worldRank, this->worldSize)));
#endif

#ifdef PROFILE_3BMDA
                // only accumulate if there are possible particle interactions to avoid div by 0
                if (std::get<1>(numParticleInteractions) > 0) {
                    hitRateDivider += std::get<1>(numParticleInteractions);
                }
#endif
            }

            if (this->worldSize > 1) {
#ifdef PROFILE_3BMDA
                // this prevents double time measurements
                std::chrono::time_point<std::chrono::system_clock> start;
                std::chrono::time_point<std::chrono::system_clock> end;
                start = std::chrono::system_clock::now();

                MPI_Barrier(this->simulation->GetTopology()->GetComm());

                end = std::chrono::system_clock::now();
                auto elapsed_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
                bool hasKey = this->times.count("idle");
                if (!hasKey) {
                    this->times["idle"] = std::make_pair(0, std::vector<int64_t>());
                }
                this->times["idle"].second.push_back(elapsed_time.count());
#endif
                shiftRight(b2);
                ++step;
            }
        }
        if (this->worldSize > 1) {
#ifdef PROFILE_3BMDA
            // this prevents double time measurements
            std::chrono::time_point<std::chrono::system_clock> start;
            std::chrono::time_point<std::chrono::system_clock> end;
            start = std::chrono::system_clock::now();

            MPI_Barrier(this->simulation->GetTopology()->GetComm());

            end = std::chrono::system_clock::now();
            auto elapsed_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
            bool hasKey = this->times.count("idle");
            if (!hasKey) {
                this->times["idle"] = std::make_pair(0, std::vector<int64_t>());
            }
            this->times["idle"].second.push_back(elapsed_time.count());
#endif
            shiftRight(b1);
        }
    }

#ifdef PROFILE_3BMDA
    if (hitRateDivider > 0) {
        this->hitrates.push_back((double)numParticleInteractionsAcc / (double)hitRateDivider);
    }
#endif

    this->SumUpParticles(this->b0, this->b1, this->b2);

    this->simulation->GetDecomposition()->SetMyParticles(this->b0);

    return std::tuple(numBufferInteractions, numParticleInteractionsAcc);
}

std::tuple<uint64_t, uint64_t> NATA::SimulationStep([[maybe_unused]] ForceType forceType) {
    throw std::runtime_error("not implemented");
}