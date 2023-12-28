#include "AUTA.hpp"

AUTA::AUTA() {}

AUTA::~AUTA() {}

void AUTA::Init(std::shared_ptr<Simulation> simulation) {
    Algorithm::Init(simulation);

    // in this algorithm we copy b0, as this buffer is also shifted around
    this->b0 = this->simulation->GetDecomposition()->GetMyParticles();
    this->ringTopology = (std::static_pointer_cast<RingTopology>(this->simulation->GetTopology()));
    this->leftNeighbor = ringTopology->GetLeftNeighbor();
    this->rightNeighbor = ringTopology->GetRightNeighbor();
    this->worldRank = ringTopology->GetWorldRank();
    this->worldSize = ringTopology->GetWorldSize();

    this->b0Owner = this->worldRank;
    this->b1Owner = this->worldRank;
    this->b2Owner = this->worldRank;
}

int AUTA::shiftRight(std::vector<Utility::Particle>& buf, int owner) {
#ifdef PROFILE_3BMDA
    std::chrono::time_point<std::chrono::system_clock> start;
    std::chrono::time_point<std::chrono::system_clock> end;
    start = std::chrono::system_clock::now();
#endif

    MPI_Status status;

    // https://moodle.rrze.uni-erlangen.de/pluginfile.php/13157/mod_resource/content/1/06_MPI_Advanced.pdf Page 12
    MPI_Sendrecv_replace(buf.data(), buf.size(), *this->mpiParticleType, this->rightNeighbor, owner, this->leftNeighbor,
                         MPI_ANY_TAG, this->ringTopology->GetComm(), &status);

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

    return status.MPI_TAG;
}

std::tuple<uint64_t, uint64_t> AUTA::calculateOneThirdOfInteractions(int thirdID) {
    std::vector<Utility::Particle>* b0Sorted = nullptr;
    std::vector<Utility::Particle>* b1Sorted = nullptr;
    std::vector<Utility::Particle>* b2Sorted = nullptr;
    int b0OwnerSorted = this->b0Owner;
    int b1OwnerSorted = this->b1Owner;
    int b2OwnerSorted = this->b2Owner;

    // sort buffers by owner
    if (this->b0Owner < this->b1Owner && this->b0Owner < this->b2Owner) {
        if (this->b1Owner < this->b2Owner) {
            b0Sorted = &(this->b0);
            b1Sorted = &(this->b1);
            b2Sorted = &(this->b2);
        } else {
            b0Sorted = &(this->b0);
            b1Sorted = &(this->b2);
            b2Sorted = &(this->b1);
            b0OwnerSorted = this->b0Owner;
            b1OwnerSorted = this->b2Owner;
            b2OwnerSorted = this->b1Owner;
        }
    } else if (this->b1Owner < this->b0Owner && this->b1Owner < this->b2Owner) {
        if (this->b0Owner < this->b2Owner) {
            b0Sorted = &(this->b1);
            b1Sorted = &(this->b0);
            b2Sorted = &(this->b2);
            b0OwnerSorted = this->b1Owner;
            b1OwnerSorted = this->b0Owner;
            b2OwnerSorted = this->b2Owner;
        } else {
            b0Sorted = &(this->b1);
            b1Sorted = &(this->b2);
            b2Sorted = &(this->b0);
            b0OwnerSorted = this->b1Owner;
            b1OwnerSorted = this->b2Owner;
            b2OwnerSorted = this->b0Owner;
        }
    } else if (this->b2Owner < this->b1Owner && this->b2Owner < this->b0Owner) {
        if (this->b0Owner < this->b1Owner) {
            b0Sorted = &(this->b2);
            b1Sorted = &(this->b0);
            b2Sorted = &(this->b1);
            b0OwnerSorted = this->b2Owner;
            b1OwnerSorted = this->b0Owner;
            b2OwnerSorted = this->b1Owner;
        } else {
            b0Sorted = &(this->b2);
            b1Sorted = &(this->b1);
            b2Sorted = &(this->b0);
            b0OwnerSorted = this->b2Owner;
            b1OwnerSorted = this->b1Owner;
            b2OwnerSorted = this->b0Owner;
        }
    }

    int start = thirdID * (b0Sorted->size() / 3);
    int numSteps = b0Sorted->size() / 3;

    // the last processor calculates the rest if #particles in b0Sorted is not divisable by 3
    if (thirdID == 2) {
        numSteps = b0Sorted->size() - 2 * numSteps;
    }

    return this->CalculateInteractions(*b0Sorted, *b1Sorted, *b2Sorted, b0OwnerSorted, b1OwnerSorted, b2OwnerSorted,
                                       start, numSteps);
}

std::vector<Utility::Particle>* AUTA::pickBuffer(int i) {
    switch (i) {
        case 0:
            return &this->b0;
            break;
        case 1:
            return &this->b1;
            break;
        case 2:
            return &this->b2;
            break;

        default:
            exit(1);
    }
}

int& AUTA::getBufOwner(int i) {
    switch (i) {
        case 0:
            return this->b0Owner;
            break;
        case 1:
            return this->b1Owner;
            break;
        case 2:
            return this->b2Owner;
            break;

        default:
            exit(1);
    }
}

void AUTA::sendBackParticles() {
#ifdef PROFILE_3BMDA
    std::chrono::time_point<std::chrono::system_clock> start;
    std::chrono::time_point<std::chrono::system_clock> end;
    start = std::chrono::system_clock::now();
#endif
    MPI_Request requestSend0, requestSend1, requestSend2;
    MPI_Request requestRecv0, requestRecv1, requestRecv2;

    bool b0Sent = false, b1Sent = false, b2Sent = false;

    if (this->b0Owner != this->worldRank) {
        MPI_Isend(this->b0.data(), this->b0.size(), *this->mpiParticleType, this->b0Owner, 0,
                  this->ringTopology->GetComm(), &requestSend0);
        b0Sent = true;
    }
    if (this->b1Owner != this->worldRank) {
        MPI_Isend(this->b1.data(), this->b1.size(), *this->mpiParticleType, this->b1Owner, 1,
                  this->ringTopology->GetComm(), &requestSend1);
        b1Sent = true;
    }
    if (this->b2Owner != this->worldRank) {
        MPI_Isend(this->b2.data(), this->b2.size(), *this->mpiParticleType, this->b2Owner, 2,
                  this->ringTopology->GetComm(), &requestSend2);
        b2Sent = true;
    }

    // all buffers have the same size
    int numRecv = b0.size();

    if (this->b0Owner != this->worldRank) {
        this->b0Tmp.resize(numRecv);

        MPI_Irecv(b0Tmp.data(), numRecv, *this->mpiParticleType, MPI_ANY_SOURCE, 0, this->ringTopology->GetComm(),
                  &requestRecv0);
    }
    if (this->b1Owner != this->worldRank) {
        this->b1Tmp.resize(numRecv);

        MPI_Irecv(b1Tmp.data(), numRecv, *this->mpiParticleType, MPI_ANY_SOURCE, 1, this->ringTopology->GetComm(),
                  &requestRecv1);
    }
    if (this->b2Owner != this->worldRank) {
        this->b2Tmp.resize(numRecv);

        MPI_Irecv(b2Tmp.data(), numRecv, *this->mpiParticleType, MPI_ANY_SOURCE, 2, this->ringTopology->GetComm(),
                  &requestRecv2);
    }

    if (b0Sent) {
        MPI_Wait(&requestSend0, MPI_STATUS_IGNORE);
        MPI_Wait(&requestRecv0, MPI_STATUS_IGNORE);
    }
    if (b1Sent) {
        MPI_Wait(&requestSend1, MPI_STATUS_IGNORE);
        MPI_Wait(&requestRecv1, MPI_STATUS_IGNORE);
    }
    if (b2Sent) {
        MPI_Wait(&requestSend2, MPI_STATUS_IGNORE);
        MPI_Wait(&requestRecv2, MPI_STATUS_IGNORE);
    }

    if (b0Sent) {
        this->b0 = this->b0Tmp;
        this->b0Tmp.clear();
    }

    if (b1Sent) {
        this->b1 = this->b1Tmp;
        this->b1Tmp.clear();
    }

    if (b2Sent) {
        this->b2 = this->b2Tmp;
        this->b2Tmp.clear();
    }

    this->b0Owner = this->worldRank;
    this->b1Owner = this->worldRank;
    this->b2Owner = this->worldRank;
#ifdef PROFILE_3BMDA
    end = std::chrono::system_clock::now();
    auto elapsed_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    bool hasKey = this->times.count("sendBackParticles");
    if (!hasKey) {
        this->times["sendBackParticles"] = std::make_pair(0, std::vector<int64_t>());
    }
    this->times["sendBackParticles"].second.push_back(elapsed_time.count());
#endif
}

std::tuple<uint64_t, uint64_t> AUTA::SimulationStep() {
    this->b0 = this->simulation->GetDecomposition()->GetMyParticles();

    this->numShifts = 0;

    // use assignment operator to copy vector
    b1 = b0;
    b2 = b0;

    this->b0Owner = this->worldRank;
    this->b1Owner = this->worldRank;
    this->b2Owner = this->worldRank;

    int i = 2;
    std::vector<Utility::Particle>* bi = pickBuffer(i);

    uint64_t numBufferInteractions = 0;
    uint64_t numParticleInteractionsAcc = 0;

#ifdef TESTS_3BMDA
    processed.clear();
#endif

#ifdef PROFILE_3BMDA
    this->hitrate = 0;
    uint64_t hitRateDivider = 0;
#endif

    for (int s = this->worldSize; s > 0; s -= 3) {
        for (int j = 0; j < s; ++j) {
            if (j != 0 || s != this->worldSize) {
#ifdef PROFILE_3BMDA
                // this prevents double time measurements
                std::chrono::time_point<std::chrono::system_clock> start;
                std::chrono::time_point<std::chrono::system_clock> end;
                start = std::chrono::system_clock::now();

                MPI_Barrier(this->simulation->GetTopology()->GetComm());

                end = std::chrono::system_clock::now();
#endif
                getBufOwner(i) = shiftRight(*bi, getBufOwner(i));

#ifdef PROFILE_3BMDA
                auto elapsed_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
                bool hasKey = this->times.count("idle");
                if (!hasKey) {
                    this->times["idle"] = std::make_pair(0, std::vector<int64_t>());
                }
                this->times["idle"].second.push_back(elapsed_time.count());
#endif
            }

#if defined(VLEVEL) && !defined(BENCHMARK_3BMDA) && !defined(TESTS_3BMDA) && VLEVEL > 0
            std::string message = "I'm proc " + std::to_string(simulation->GetTopology()->GetWorldRank()) +
                                  " and going to calculate interactions between buffers from proc (" +
                                  std::to_string(this->b0Owner) + ", " + std::to_string(this->b1Owner) + ", " +
                                  std::to_string(this->b2Owner) + ")";
            MPIReporter::instance()->StoreMessage(this->simulation->GetTopology()->GetWorldRank(), message);
#endif

            std::tuple<uint64_t, uint64_t> numParticleInteractions =
                this->CalculateInteractions(this->b0, this->b1, this->b2, this->b0Owner, this->b1Owner, this->b2Owner);
            numParticleInteractionsAcc += std::get<0>(numParticleInteractions);

            numBufferInteractions++;
#ifdef TESTS_3BMDA
            processed.push_back(Utility::Triplet(this->b0Owner, this->b1Owner, this->b2Owner));
#endif
#ifdef PROFILE_3BMDA
            // only accumulate if there are possible particle interactions to avoid div by 0
            if (std::get<1>(numParticleInteractions) > 0) {
                hitRateDivider += std::get<1>(numParticleInteractions);
            }
#endif
        }
        i = (i + 1) % 3;
        bi = pickBuffer(i);
    }
    if (this->worldSize % 3 == 0) {
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
        getBufOwner(i) = shiftRight(*bi, getBufOwner(i));

        int thirdID = this->worldRank / (this->worldSize / 3);

        // Calculate one third of the interactions

#if defined(VLEVEL) && !defined(BENCHMARK_3BMDA) && !defined(TESTS_3BMDA) && VLEVEL > 0
        std::string message = "I'm proc " + std::to_string(simulation->GetTopology()->GetWorldRank()) +
                              " and going to calculate interactions between buffers from proc (" +
                              std::to_string(this->b0Owner) + ", " + std::to_string(this->b1Owner) + ", " +
                              std::to_string(this->b2Owner) + ")";
        MPIReporter::instance()->StoreMessage(this->simulation->GetTopology()->GetWorldRank(), message);
#endif

        numBufferInteractions++;
        std::tuple<uint64_t, uint64_t> numParticleInteractions = calculateOneThirdOfInteractions(thirdID);
        numParticleInteractionsAcc += std::get<0>(numParticleInteractions);

#ifdef TESTS_3BMDA
        // TESTS_3BMDA is defined
        processed.push_back(Utility::Triplet(getBufOwner(0), getBufOwner(1), getBufOwner(2)));
#endif

#ifdef PROFILE_3BMDA
        // only accumulate if there are possible particle interactions to avoid div by 0
        if (std::get<1>(numParticleInteractions) > 0) {
            hitRateDivider += std::get<1>(numParticleInteractions);
        }
#endif
    }

    // send back to owner
    if (this->worldSize > 1) {
        sendBackParticles();
    }

#ifdef PROFILE_3BMDA
    if (hitRateDivider > 0) {
        this->hitrates.push_back((double)numParticleInteractionsAcc / (double)hitRateDivider);
    }
#endif

    // sum up particles
    this->SumUpParticles(this->b0, this->b1, this->b2);

    this->simulation->GetDecomposition()->SetMyParticles(this->b0);

    return std::tuple(numBufferInteractions, numParticleInteractionsAcc);
}

std::tuple<uint64_t, uint64_t> AUTA::SimulationStep([[maybe_unused]] ForceType forceType) {
    throw std::runtime_error("not implemented");
}