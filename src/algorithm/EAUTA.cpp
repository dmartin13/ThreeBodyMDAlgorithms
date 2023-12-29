#include "EAUTA.hpp"

EAUTA::EAUTA() {}

EAUTA::~EAUTA() {}

void EAUTA::Init(std::shared_ptr<Simulation> simulation) {
    Algorithm::Init(simulation);

    // in this algorithm we copy b0, as this buffer is also shifted around
    this->b1 = this->simulation->GetDecomposition()->GetMyParticles();
    this->ringTopology = (std::static_pointer_cast<RingTopology>(this->simulation->GetTopology()));
    this->leftNeighbor = ringTopology->GetLeftNeighbor();
    this->rightNeighbor = ringTopology->GetRightNeighbor();
    this->worldRank = ringTopology->GetWorldRank();
    this->worldSize = ringTopology->GetWorldSize();

    this->b0Owner = this->worldRank;
    this->b1Owner = this->worldRank;
    this->b2Owner = this->worldRank;
}

int EAUTA::shiftRight(std::vector<Utility::Particle>& buf, int owner) {
    MPI_Status status;

    MPI_Sendrecv_replace(buf.data(), buf.size(), *this->mpiParticleType, this->rightNeighbor, owner, this->leftNeighbor,
                         MPI_ANY_TAG, this->ringTopology->GetComm(), &status);

    ++numShifts;

    // returns the owner
    return status.MPI_TAG;
}

int EAUTA::shiftLeft(std::vector<Utility::Particle>& buf, int owner) {
    MPI_Status status;

    MPI_Sendrecv_replace(buf.data(), buf.size(), *this->mpiParticleType, this->leftNeighbor, owner, this->rightNeighbor,
                         MPI_ANY_TAG, this->ringTopology->GetComm(), &status);

    ++numShifts;

    // returns the owner
    return status.MPI_TAG;
}

std::tuple<uint64_t, uint64_t> EAUTA::calculateOneThirdOfInteractions(int thirdID) {
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

std::tuple<uint64_t, uint64_t> EAUTA::calculateOneHalfOfInteractions(std::vector<Utility::Particle>& b0,
                                                                     std::vector<Utility::Particle>& b1, int b0Owner,
                                                                     int b1Owner, int halfID) {
    std::vector<Utility::Particle>* b0Sorted = nullptr;
    std::vector<Utility::Particle>* b1Sorted = nullptr;
    int b0OwnerSorted;
    int b1OwnerSorted;

    // sort buffers by owner
    if (b0Owner < b1Owner) {
        b0Sorted = &b0;
        b1Sorted = &b1;
        b0OwnerSorted = b0Owner;
        b1OwnerSorted = b1Owner;
    } else {
        b0Sorted = &b1;
        b1Sorted = &b0;
        b0OwnerSorted = b1Owner;
        b1OwnerSorted = b0Owner;
    }

    int start = halfID * (b0Sorted->size() / 2);
    int numSteps = b0Sorted->size() / 2;

    // the last processor calculates the rest if #particles in b0Sorted is not divisable by 2
    if (halfID == 1) {
        numSteps = b0Sorted->size() - numSteps;
    }

    return this->calculatePairwiseInteractions(*b0Sorted, *b1Sorted, b0OwnerSorted, b1OwnerSorted, start, numSteps, 0,
                                               Eigen::Array3d{0, 0, 0});
}

std::vector<Utility::Particle>* EAUTA::pickBuffer(int i) {
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

int& EAUTA::getBufOwner(int i) {
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

void EAUTA::sendBackParticles() {
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

    if (this->b0Owner != this->worldRank) {
        int numRecv = b0.size();
        this->b0Tmp.resize(numRecv);
        MPI_Irecv(b0Tmp.data(), numRecv, *this->mpiParticleType, MPI_ANY_SOURCE, 0, this->ringTopology->GetComm(),
                  &requestRecv0);
    }
    if (this->b1Owner != this->worldRank) {
        int numRecv = b1.size();
        this->b1Tmp.resize(numRecv);
        MPI_Irecv(b1Tmp.data(), numRecv, *this->mpiParticleType, MPI_ANY_SOURCE, 1, this->ringTopology->GetComm(),
                  &requestRecv1);
    }
    if (this->b2Owner != this->worldRank) {
        int numRecv = b2.size();
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
}

std::tuple<uint64_t, uint64_t> EAUTA::SimulationStep() { return SimulationStep(ForceType::TwoAndThreeBody); }

std::tuple<uint64_t, uint64_t> EAUTA::SimulationStep(ForceType forceType) {
    // TODO: avoid this copy
    this->b1 = this->simulation->GetDecomposition()->GetMyParticles();

    this->b1Owner = this->worldRank;

    this->numShifts = 0;

    uint64_t numBufferInteractions = 0;
    uint64_t numParticleInteractionsAcc = 0;

    // special case for 1 processor
    if (worldSize == 1) {
        if (forceType == ForceType::ThreeBody or forceType == ForceType::TwoAndThreeBody) {
            CalculateInteractions(this->b1, this->b1, this->b1, this->b1Owner, this->b1Owner, this->b1Owner);
        }
        if (forceType == ForceType::TwoBody or forceType == ForceType::TwoAndThreeBody) {
            CalculatePairwiseInteractions(b1, b1, b1Owner, b1Owner);
        }
        this->simulation->GetDecomposition()->SetMyParticles(this->b1);
        return std::tuple(numBufferInteractions, numParticleInteractionsAcc);
    }

    // special case for 2 processors
    if (worldSize == 2) {
        b2 = b1;

        this->b2Owner = shiftLeft(b2, worldRank);

        if (forceType == ForceType::ThreeBody or forceType == ForceType::TwoAndThreeBody) {
            CalculateInteractions(this->b1, this->b1, this->b1, this->b1Owner, this->b1Owner, this->b1Owner);
            CalculateInteractions(this->b1, this->b1, this->b2, this->b1Owner, this->b1Owner, this->b2Owner);
        }
        if (forceType == ForceType::TwoBody or forceType == ForceType::TwoAndThreeBody) {
            CalculatePairwiseInteractions(b1, b1, b1Owner, b1Owner);

            // Calculate one half of interactions
            calculateOneHalfOfInteractions(b1, b2, b1Owner, b2Owner, worldRank);
        }

        sendBackParticles();

        SumUpParticles(this->b1, this->b2);

        simulation->GetDecomposition()->SetMyParticles(this->b1);

        return std::tuple(numBufferInteractions, numParticleInteractionsAcc);
    }

    // special case for 3 processors
    if (worldSize == 3) {
        b0 = b1;
        b2 = b1;

        this->b0Owner = shiftRight(b0, worldRank);
        this->b2Owner = shiftLeft(b2, worldRank);

        if (forceType == ForceType::ThreeBody or forceType == ForceType::TwoAndThreeBody) {
            CalculateInteractions(this->b1, this->b1, this->b1, this->b1Owner, this->b1Owner, this->b1Owner);
            CalculateInteractions(this->b1, this->b1, this->b2, this->b1Owner, this->b1Owner, this->b2Owner);
            CalculateInteractions(this->b0, this->b0, this->b2, this->b0Owner, this->b0Owner, this->b2Owner);

            // special case
            int thirdID = this->worldRank / (this->worldSize / 3);
            // Calculate one third of the interactions between b0, b1, b2
            calculateOneThirdOfInteractions(thirdID);
        }

        if (forceType == ForceType::TwoBody or forceType == ForceType::TwoAndThreeBody) {
            CalculatePairwiseInteractions(b1, b1, b1Owner, b1Owner);
            CalculatePairwiseInteractions(b1, b2, b1Owner, b2Owner);
        }

        sendBackParticles();

        SumUpParticles(this->b0, this->b1, this->b2);
        simulation->GetDecomposition()->SetMyParticles(this->b0);

        return std::tuple(numBufferInteractions, numParticleInteractionsAcc);
    }

    // use assignment operator to copy vector
    b0 = b1;
    b2 = b1;

    this->b0Owner = shiftRight(b0, worldRank);
    this->b2Owner = shiftLeft(b2, worldRank);

    int i = 0;
    std::vector<Utility::Particle>* bi = pickBuffer(i);

    // The two-body interactions are integrated in sequence into the first substeps of the first phase in the AUTA
    // algorithm. Since we have (worldSize + 2 - 1) choose 2 unique pairs distributed over worldSize processors, we can
    // stop two-body calculation after we have done the required number of steps.
    size_t numStepsPairwise = Utility::BinomialCoefficient(worldSize + 2 - 1, 2) / worldSize;
    size_t pairCounter = 0;

#ifdef TESTS_3BMDA
    processed.clear();
#endif

    for (int s = this->worldSize - 3; s > 0; s -= 3) {
        for (int j = 0; j < s; ++j) {
            if (j != 0 || s != this->worldSize - 3) {
                getBufOwner(i) = shiftRight(*bi, getBufOwner(i));
            } else {
                // in the very first step (substep 1 of phase 1) we incorporate some 3-body and 2-body interactions we
                // can calculate from the buffers we already have
                if (forceType == ForceType::ThreeBody or forceType == ForceType::TwoAndThreeBody) {
                    CalculateInteractions(this->b1, this->b1, this->b1, this->b1Owner, this->b1Owner, this->b1Owner);
                    CalculateInteractions(this->b1, this->b1, this->b2, this->b1Owner, this->b1Owner, this->b2Owner);
                    CalculateInteractions(this->b0, this->b0, this->b2, this->b0Owner, this->b0Owner, this->b2Owner);
                }

                if (forceType == ForceType::TwoBody or forceType == ForceType::TwoAndThreeBody) {
                    // MPIReporter::instance()->StoreMessage(
                    //     worldRank, "calculates (" + std::to_string(b1Owner) + ", " + std::to_string(b1Owner) + "), ("
                    //     +
                    //                    std::to_string(b0Owner) + ", " + std::to_string(b1Owner) + ")");
                    CalculatePairwiseInteractions(b1, b1, b1Owner, b1Owner);
                    CalculatePairwiseInteractions(b0, b1, b0Owner, b1Owner);
                    pairCounter += 2;
                }
            }
            if (s == this->worldSize - 3) {
                // See Embedded AUT Algorithm paper: In the first phase we incorporate also some additional 3-body
                // calculations after buffer 0 is shifted
                if (forceType == ForceType::ThreeBody or forceType == ForceType::TwoAndThreeBody) {
                    CalculateInteractions(this->b0, this->b1, this->b1, this->b0Owner, this->b1Owner, this->b1Owner);
                }
            }

            // normal interactions between all three buffers
            if (forceType == ForceType::ThreeBody or forceType == ForceType::TwoAndThreeBody) {
                CalculateInteractions(this->b0, this->b1, this->b2, this->b0Owner, this->b1Owner, this->b2Owner);
            }

            // calculate pair interactions until we have all necessary pairs (note: if the number of processors is
            // divisible by 2, we need one extra round where interactions are distributed over processors)
            if (forceType == ForceType::TwoBody or forceType == ForceType::TwoAndThreeBody) {
                if ((pairCounter < numStepsPairwise) or (worldSize % 2 == 0 and pairCounter == numStepsPairwise)) {
                    if (pairCounter < numStepsPairwise) {
                        // MPIReporter::instance()->StoreMessage(worldRank, "calculates no special case (" +
                        //                                                      std::to_string(b0Owner) + ", " +
                        //                                                      std::to_string(b2Owner) + ")");

                        CalculatePairwiseInteractions(b0, b2, b0Owner, b2Owner);
                    } else {
                        // MPIReporter::instance()->StoreMessage(
                        //     worldRank,
                        //     "calculates special case (" + std::to_string(b0Owner) + ", " + std::to_string(b2Owner) +
                        //     ")");
                        // special case for pairs
                        int halfID = this->worldRank / (this->worldSize / 2);
                        calculateOneHalfOfInteractions(b0, b2, b0Owner, b2Owner, halfID);
                    }
                    ++pairCounter;
                }
            }

            // since two-body interactions are fully embedded into the first phase of the algorithm we can stop the
            // algorithm here, if we have a non-respa iteration
            if (forceType == ForceType::TwoBody and
                ((worldSize % 2) == 0 ? pairCounter > numStepsPairwise : pairCounter >= numStepsPairwise)) {
                // send back to owner
                sendBackParticles();
                // sum up particles
                this->SumUpParticles(this->b0, this->b1, this->b2);
                this->simulation->GetDecomposition()->SetMyParticles(this->b0);
                return std::tuple(numBufferInteractions, numParticleInteractionsAcc);
            }

#ifdef TESTS_3BMDA
            processed.push_back(Utility::Triplet(this->b0Owner, this->b1Owner, this->b2Owner));
#endif
        }
        i = (i + 1) % 3;
        bi = pickBuffer(i);
    }
    if (this->worldSize % 3 == 0 and (forceType == ForceType::ThreeBody or forceType == ForceType::TwoAndThreeBody)) {
        getBufOwner(i) = shiftRight(*bi, getBufOwner(i));

        int thirdID = this->worldRank / (this->worldSize / 3);

        // Calculate one third of the interactions between b0, b1, b2
        calculateOneThirdOfInteractions(thirdID);

#ifdef TESTS_3BMDA
        // TESTS_3BMDA is defined
        processed.push_back(Utility::Triplet(getBufOwner(0), getBufOwner(1), getBufOwner(2)));
#endif
    }

    // send back to owner
    sendBackParticles();

    // sum up particles
    this->SumUpParticles(this->b0, this->b1, this->b2);

    this->simulation->GetDecomposition()->SetMyParticles(this->b0);

    return std::tuple(numBufferInteractions, numParticleInteractionsAcc);
}