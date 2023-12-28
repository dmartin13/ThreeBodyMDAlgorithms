#include "P3BCA.hpp"

P3BCA::P3BCA(double cutoff) : cutoff(cutoff) {}

P3BCA::~P3BCA() {}

void P3BCA::Init(std::shared_ptr<Simulation> simulation) {
    Algorithm::Init(simulation);

    this->b0 = this->simulation->GetDecomposition()->GetMyParticles();
    this->cartTopology = std::static_pointer_cast<CartTopology>(this->simulation->GetTopology());

    std::shared_ptr<RegularGridDecomposition> decomposition =
        std::static_pointer_cast<RegularGridDecomposition>(this->simulation->GetDecomposition());

    this->physicalDomainSize = decomposition->GetPhysicalDomainSize();

    this->dimX = decomposition->GetDimX();
    this->dimY = decomposition->GetDimY();
    this->dimZ = decomposition->GetDimZ();

    Eigen::Array3d localCellSize = decomposition->GetCellSize();
    double cellSizeX = localCellSize[0];
    double cellSizeY = localCellSize[1];
    double cellSizeZ = localCellSize[2];

    this->numDims = this->cartTopology->GetCartRank().GetDimensions();

    // all boxed are included that are fully or partially covered by the cutoff distance
    if (cartTopology->GetWorldSize() == 1) {
        this->nCbX = 0;
    } else if (cartTopology->GetWorldSize() == 2) {
        this->nCbX = 1;
    } else {
        this->nCbX = (int)std::ceil((this->cutoff / cellSizeX));
    }
    if (numDims > 1) {
        this->nCbY = (int)std::ceil((this->cutoff / cellSizeY));
    } else {
        this->nCbY = 0;
    }
    if (numDims > 2) {
        this->nCbZ = (int)std::ceil((this->cutoff / cellSizeZ));
    } else {
        this->nCbZ = 0;
    }

    this->worldRank = this->cartTopology->GetWorldRank();

    calcSteps(numDims);

    if (this->worldRank == 0) {
        if (this->cartTopology->GetWorldSize() > 2 &&
            (this->nCbX < 1 || (double)this->nCbX >= ((double)this->dimX / 2.0))) {
            std::cerr << "invalid cutoff: " << this->cutoff << ", with " << this->nCbX
                      << " cutoff boxes, dimX: " << this->dimX << ", and a cell size of " << cellSizeX << std::endl;
            exit(1);
        }
        if (this->numDims > 1 && (this->nCbY < 1 || (double)this->nCbY >= ((double)this->dimY / 2.0))) {
            std::cerr << "invalid cutoff: " << this->cutoff << ", with " << this->nCbY
                      << " cutoff boxes, dimY: " << this->dimY << ", and a cell size of " << cellSizeY << std::endl;
            exit(1);
        }
        if (this->numDims > 2 && (this->nCbZ < 1 || (double)this->nCbZ >= ((double)this->dimZ / 2.0))) {
            std::cerr << "invalid cutoff: " << this->cutoff << ", with " << this->nCbZ
                      << " cutoff boxes, dimZ: " << this->dimZ << ", and a cell size of " << cellSizeZ << std::endl;
            exit(1);
        }

#if defined(VLEVEL) && !defined(BENCHMARK_3BMDA) && !defined(TESTS_3BMDA)
        std::string message1 = "numDims: " + std::to_string(numDims);

        std::string message2 = "specified cutoff: " + std::to_string(this->cutoff) + ", calculated cutoff boxes: [" +
                               std::to_string(this->nCbX) + ", " + std::to_string(this->nCbY) + ", " +
                               std::to_string(this->nCbZ) + "], cell size: [" + std::to_string(cellSizeX) + ", " +
                               std::to_string(cellSizeY) + ", " + std::to_string(cellSizeZ) + "]";

        MPIReporter::instance()->StoreMessage(this->simulation->GetTopology()->GetWorldRank(), message1);
        MPIReporter::instance()->StoreMessage(this->simulation->GetTopology()->GetWorldRank(), message2);
#endif
    }
}

void P3BCA::calcSteps(int dimension) {
    switch (dimension) {
        // do just one simulation step that calculates all interactions if we have only one processor
        case 1:
            this->numSteps = this->nCbX + 1;
            break;
        case 2:
            this->numSteps = 2 * this->nCbX * this->nCbY + this->nCbX + this->nCbY + 1;
            break;
        case 3:
            this->numSteps = 1 + this->nCbZ + this->nCbY * (this->nCbZ * 2 + 1) +
                             this->nCbX * ((this->nCbY * 2 + 1) * (this->nCbZ * 2 + 1));
            break;

        default:
            break;
    }
}

void P3BCA::schedule1D(int i, int& myCartRank, int& src) { src = Utility::mod(i + myCartRank, this->dimX); }

void P3BCA::calcDestFromSrc1D(int& myCartRank, int& src, int& dst) {
    int shiftedSrc;
    shiftedSrc = Utility::mod(src - myCartRank, this->dimX);

    dst = Utility::mod(myCartRank - shiftedSrc, this->dimX);
}

void P3BCA::schedule1DHelper(int i3, int& cartRank, int& src, int& dst, int& diff) {
    int initI3 = i3;

    src = Utility::mod((i3 + 1) + cartRank, this->dimX);

    calcDestFromSrc1D(cartRank, src, dst);
    calcDiff1D(cartRank, src, diff, initI3 + 1);
}

void P3BCA::calcDiff1D(int& cartRank, int& src, int& diff, int i) {
    int oldSrc;
    int iOld = i - 1;

    oldSrc = Utility::mod(iOld + cartRank, this->dimX);

    diff = oldSrc - src;
}

void P3BCA::schedule2D(int i, std::array<int, 2>& myCartRank, std::array<int, 2>& src) {
    int cutoffLenY = 2 * this->nCbY + 1;

    int g_y_i = ((i + this->nCbY) % cutoffLenY) - this->nCbY;
    int g_x_i = (i + this->nCbY) / cutoffLenY;

    src[0] = Utility::mod(g_x_i + myCartRank[0], this->dimX);
    src[1] = Utility::mod(g_y_i + myCartRank[1], this->dimY);
}

void P3BCA::calcDestFromSrc2D(std::array<int, 2>& myCartRank, std::array<int, 2>& src, std::array<int, 2>& dst) {
    std::array<int, 2> shiftedSrc;
    shiftedSrc[0] = Utility::mod(src[0] - myCartRank[0], this->dimX);
    shiftedSrc[1] = Utility::mod(src[1] - myCartRank[1], this->dimY);

    dst[0] = Utility::mod(myCartRank[0] - shiftedSrc[0], this->dimX);
    dst[1] = Utility::mod(myCartRank[1] - shiftedSrc[1], this->dimY);
}

void P3BCA::schedule2DHelper(int i2, int& i3, std::array<int, 2>& cartRank, std::array<int, 2>& src,
                             std::array<int, 2>& dst, std::array<int, 2>& diff) {
    int myCoords[2];
    MPI_Cart_coords(this->cartTopology->GetComm(), this->worldRank, 2, myCoords);

    std::array<int, 2> b1InitOwner;
    std::array<int, 2> srcTmp;
    std::array<int, 2> srcTmp2;

    int initI3 = i3;

    do {
        // find the rank of the processor that initially owned the paticles which are in our b1 at the moment
        schedule2D(i2, cartRank, b1InitOwner);

        // from our processors rank get the processor if we would move to the next processor after i3 steps
        schedule2D(i3 + 1, cartRank, srcTmp);

        // from b1InitOwner move to the next processor
        schedule2D(i3 + 1 - i2, b1InitOwner, srcTmp2);

    } while (srcTmp != srcTmp2 && i3 < this->numSteps && i3++);

    src = srcTmp;

    calcDestFromSrc2D(cartRank, src, dst);

    calcDiff2D(cartRank, src, diff, initI3 + 1);
}

// a mapping from an int to a cartesian coordinate (int, int, int)
void P3BCA::schedule3D(int i, std::array<int, 3>& myCartRank, std::array<int, 3>& src) {
    int cutoffLenY = 2 * this->nCbY + 1;
    int cutoffLenZ = 2 * this->nCbZ + 1;

    int g_z_i = ((i + this->nCbZ) % cutoffLenZ) - this->nCbZ;
    int g_y_i = ((((i + this->nCbZ) / cutoffLenZ) + this->nCbY) % cutoffLenY) - this->nCbY;
    int g_x_i = (i + ((cutoffLenZ * cutoffLenY) / 2)) / (cutoffLenZ * cutoffLenY);

    src[0] = Utility::mod(g_x_i + myCartRank[0], this->dimX);
    src[1] = Utility::mod(g_y_i + myCartRank[1], this->dimY);
    src[2] = Utility::mod(g_z_i + myCartRank[2], this->dimZ);
}

void P3BCA::calcDestFromSrc3D(std::array<int, 3>& myCartRank, std::array<int, 3>& src, std::array<int, 3>& dst) {
    std::array<int, 3> shiftedSrc;
    shiftedSrc[0] = Utility::mod(src[0] - myCartRank[0], this->dimX);
    shiftedSrc[1] = Utility::mod(src[1] - myCartRank[1], this->dimY);
    shiftedSrc[2] = Utility::mod(src[2] - myCartRank[2], this->dimZ);

    dst[0] = Utility::mod(myCartRank[0] - shiftedSrc[0], this->dimX);
    dst[1] = Utility::mod(myCartRank[1] - shiftedSrc[1], this->dimY);
    dst[2] = Utility::mod(myCartRank[2] - shiftedSrc[2], this->dimZ);
}

void P3BCA::calcDiff3D(std::array<int, 3>& cartRank, std::array<int, 3>& src, std::array<int, 3>& diff, int i) {
    int cutoffLenY = 2 * this->nCbY + 1;
    int cutoffLenZ = 2 * this->nCbZ + 1;

    std::array<int, 3> oldSrc;
    int iOld = i - 1;

    int g_z_iOld = ((iOld + this->nCbZ) % cutoffLenZ) - this->nCbZ;
    int g_y_iOld = ((((iOld + this->nCbZ) / cutoffLenZ) + this->nCbY) % cutoffLenY) - this->nCbY;
    int g_x_iOld = (iOld + ((cutoffLenZ * cutoffLenY) / 2)) / (cutoffLenZ * cutoffLenY);

    oldSrc[0] = Utility::mod(g_x_iOld + cartRank[0], this->dimX);
    oldSrc[1] = Utility::mod(g_y_iOld + cartRank[1], this->dimY);
    oldSrc[2] = Utility::mod(g_z_iOld + cartRank[2], this->dimZ);

    diff[0] = oldSrc[0] - src[0];
    diff[1] = oldSrc[1] - src[1];
    diff[2] = oldSrc[2] - src[2];
}

void P3BCA::calcDiff2D(std::array<int, 2>& cartRank, std::array<int, 2>& src, std::array<int, 2>& diff, int i) {
    int cutoffLenY = 2 * this->nCbY + 1;

    std::array<int, 2> oldSrc;
    int iOld = i - 1;

    int g_y_iOld = ((iOld + this->nCbY) % cutoffLenY) - this->nCbY;
    int g_x_iOld = (iOld + this->nCbY) / cutoffLenY;

    oldSrc[0] = Utility::mod(g_x_iOld + cartRank[0], this->dimX);
    oldSrc[1] = Utility::mod(g_y_iOld + cartRank[1], this->dimY);

    diff[0] = oldSrc[0] - src[0];
    diff[1] = oldSrc[1] - src[1];
}

void P3BCA::schedule3DHelper(int i2, int& i3, std::array<int, 3>& cartRank, std::array<int, 3>& src,
                             std::array<int, 3>& dst, std::array<int, 3>& diff) {
    int myCoords[3];
    MPI_Cart_coords(this->cartTopology->GetComm(), this->worldRank, 3, myCoords);

    std::array<int, 3> b1InitOwner;
    std::array<int, 3> srcTmp;
    std::array<int, 3> srcTmp2;

    int initI3 = i3;

    do {
        // find the rank of the processor that initially owned the paticles which are in our b1 at the moment
        schedule3D(i2, cartRank, b1InitOwner);

        // from our processors rank get the processor if we would move to the next processor after i3 steps
        schedule3D(i3 + 1, cartRank, srcTmp);

        // from b1InitOwner move to the next processor
        schedule3D(i3 + 1 - i2, b1InitOwner, srcTmp2);

    } while (srcTmp != srcTmp2 && i3 < this->numSteps && i3++);

    src = srcTmp;

    calcDestFromSrc3D(cartRank, src, dst);

    calcDiff3D(cartRank, src, diff, initI3 + 1);
}

void P3BCA::handleOffsetVector3D(std::array<int, 3>& nextSrcRank, std::array<int, 3>& nextDstRank,
                                 std::array<int, 3>& offsetVector, std::array<int, 3>& diff,
                                 std::array<int, 3>& coordsSrc, std::array<int, 3>& coordsDst) {
    coordsSrc[0] = Utility::mod(nextSrcRank[0] + offsetVector[0], this->dimX);
    coordsSrc[1] = Utility::mod(nextSrcRank[1] + offsetVector[1], this->dimY);
    coordsSrc[2] = Utility::mod(nextSrcRank[2] + offsetVector[2], this->dimZ);

    coordsDst[0] = Utility::mod(nextDstRank[0] + (-offsetVector[0]), this->dimX);
    coordsDst[1] = Utility::mod(nextDstRank[1] + (-offsetVector[1]), this->dimY);
    coordsDst[2] = Utility::mod(nextDstRank[2] + (-offsetVector[2]), this->dimZ);

    // adjust the offset vector for the next iteration
    offsetVector[0] += diff[0];
    offsetVector[1] += diff[1];
    offsetVector[2] += diff[2];
}

void P3BCA::handleOffsetVector2D(std::array<int, 2>& nextSrcRank, std::array<int, 2>& nextDstRank,
                                 std::array<int, 2>& offsetVector, std::array<int, 2>& diff,
                                 std::array<int, 2>& coordsSrc, std::array<int, 2>& coordsDst) {
    coordsSrc[0] = Utility::mod(nextSrcRank[0] + offsetVector[0], this->dimX);
    coordsSrc[1] = Utility::mod(nextSrcRank[1] + offsetVector[1], this->dimY);

    coordsDst[0] = Utility::mod(nextDstRank[0] + (-offsetVector[0]), this->dimX);
    coordsDst[1] = Utility::mod(nextDstRank[1] + (-offsetVector[1]), this->dimY);

    // adjust the offset vector for the next iteration
    offsetVector[0] += diff[0];
    offsetVector[1] += diff[1];
}

void P3BCA::handleOffsetVector1D(int& nextSrcRank, int& nextDstRank, int& offsetVector, int& diff, int& coordsSrc,
                                 int& coordsDst) {
    coordsSrc = Utility::mod(nextSrcRank + offsetVector, this->dimX);
    coordsDst = Utility::mod(nextDstRank + (-offsetVector), this->dimX);

    // adjust the offset vector for the next iteration
    offsetVector += diff;
}

int P3BCA::shiftLeft(std::vector<Utility::Particle>& buf, int owner, int& nextSrcRank, int& nextDstRank,
                     int& offsetVector, int& diff) {
    int coordsSrc;
    int coordsDst;

    handleOffsetVector1D(nextSrcRank, nextDstRank, offsetVector, diff, coordsSrc, coordsDst);

    int src;
    int dest;
    MPI_Cart_rank(this->cartTopology->GetComm(), &coordsSrc, &src);
    MPI_Cart_rank(this->cartTopology->GetComm(), &coordsDst, &dest);

    return mpiShift(buf, owner, src, dest);
}

int P3BCA::shiftLeft(std::vector<Utility::Particle>& buf, int owner, std::array<int, 2>& nextSrcRank,
                     std::array<int, 2>& nextDstRank, std::array<int, 2>& offsetVector, std::array<int, 2>& diff) {
    std::array<int, 2> coordsSrc;
    std::array<int, 2> coordsDst;

    handleOffsetVector2D(nextSrcRank, nextDstRank, offsetVector, diff, coordsSrc, coordsDst);

    int src;
    int dest;
    MPI_Cart_rank(this->cartTopology->GetComm(), coordsSrc.data(), &src);
    MPI_Cart_rank(this->cartTopology->GetComm(), coordsDst.data(), &dest);

    return mpiShift(buf, owner, src, dest);
}

int P3BCA::shiftLeft(std::vector<Utility::Particle>& buf, int owner, std::array<int, 3>& nextSrcRank,
                     std::array<int, 3>& nextDstRank, std::array<int, 3>& offsetVector, std::array<int, 3>& diff) {
    std::array<int, 3> coordsSrc;
    std::array<int, 3> coordsDst;

    handleOffsetVector3D(nextSrcRank, nextDstRank, offsetVector, diff, coordsSrc, coordsDst);

    int src;
    int dest;
    MPI_Cart_rank(this->cartTopology->GetComm(), coordsSrc.data(), &src);
    MPI_Cart_rank(this->cartTopology->GetComm(), coordsDst.data(), &dest);

    return mpiShift(buf, owner, src, dest);
}

int P3BCA::mpiShift(std::vector<Utility::Particle>& buf, int owner, int src, int dst) {
#ifdef PROFILE_3BMDA
    std::chrono::time_point<std::chrono::system_clock> start;
    std::chrono::time_point<std::chrono::system_clock> end;
    start = std::chrono::system_clock::now();
#endif

    this->tmpRecv.clear();

    int numRecv;
    MPI_Status status;
    MPI_Request requestSend;
    MPI_Request requestRecv;

    MPI_Isend(buf.data(), buf.size(), *this->mpiParticleType, dst, owner, this->cartTopology->GetComm(), &requestSend);

    MPI_Probe(src, MPI_ANY_TAG, this->cartTopology->GetComm(), &status);
    MPI_Get_count(&status, *this->simulation->GetMPIParticleType(), &numRecv);

    this->tmpRecv.resize(numRecv);

    MPI_Irecv(this->tmpRecv.data(), numRecv, *this->mpiParticleType, src, status.MPI_TAG, this->cartTopology->GetComm(),
              &requestRecv);

    MPI_Wait(&requestSend, MPI_STATUS_IGNORE);
    MPI_Wait(&requestRecv, MPI_STATUS_IGNORE);

    // assign particles to buf
    buf = this->tmpRecv;

#ifdef PROFILE_3BMDA
    end = std::chrono::system_clock::now();
    auto elapsed_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    bool hasKey = this->times.count("mpiShift");
    if (!hasKey) {
        this->times["mpiShift"] = std::make_pair(0, std::vector<int64_t>());
    }
    this->times["mpiShift"].second.push_back(elapsed_time.count());
#endif

    this->numShifts++;

    return status.MPI_TAG;
}

void P3BCA::sendBackParticles() {
#ifdef PROFILE_3BMDA
    std::chrono::time_point<std::chrono::system_clock> start;
    std::chrono::time_point<std::chrono::system_clock> end;
    start = std::chrono::system_clock::now();
#endif
    MPI_Request requestSend1, requestSend2;
    MPI_Request requestRecv1, requestRecv2;
    MPI_Status statusRecv1, statusRecv2;
    bool b1Sent = false, b2Sent = false;

    if (this->b1Owner != this->worldRank) {
        MPI_Isend(this->b1.data(), this->b1.size(), *this->mpiParticleType, this->b1Owner, 1,
                  this->cartTopology->GetComm(), &requestSend1);
        b1Sent = true;
    }
    if (this->b2Owner != this->worldRank) {
        MPI_Isend(this->b2.data(), this->b2.size(), *this->mpiParticleType, this->b2Owner, 2,
                  this->cartTopology->GetComm(), &requestSend2);
        b2Sent = true;
    }

#if defined(VLEVEL) && !defined(BENCHMARK_3BMDA) && !defined(TESTS_3BMDA) && VLEVEL > 0
    if (this->b1Owner != this->worldRank) {
        std::string message = "I'm proc " + std::to_string(simulation->GetTopology()->GetWorldRank()) +
                              " and going to send back particles in b1 to " + std::to_string(this->b1Owner);
        MPIReporter::instance()->StoreMessage(this->simulation->GetTopology()->GetWorldRank(), message);
    }

    if (this->b2Owner != this->worldRank) {
        std::string message = "I'm proc " + std::to_string(simulation->GetTopology()->GetWorldRank()) +
                              " and going to send back particles in b2 to " + std::to_string(this->b2Owner);
        MPIReporter::instance()->StoreMessage(this->simulation->GetTopology()->GetWorldRank(), message);
    }

#endif

    int numRecv1;
    int numRecv2;

    if (this->b1Owner != this->worldRank) {
        MPI_Probe(MPI_ANY_SOURCE, 1, this->cartTopology->GetComm(), &statusRecv1);

        MPI_Get_count(&statusRecv1, *this->simulation->GetMPIParticleType(), &numRecv1);

        this->b1Tmp.resize(numRecv1);

        MPI_Irecv(b1Tmp.data(), numRecv1, *this->mpiParticleType, statusRecv1.MPI_SOURCE, 1,
                  this->cartTopology->GetComm(), &requestRecv1);
    }
    if (this->b2Owner != this->worldRank) {
        MPI_Probe(MPI_ANY_SOURCE, 2, this->cartTopology->GetComm(), &statusRecv2);

        MPI_Get_count(&statusRecv2, *this->simulation->GetMPIParticleType(), &numRecv2);

        this->b2Tmp.resize(numRecv2);

        MPI_Irecv(b2Tmp.data(), numRecv2, *this->mpiParticleType, statusRecv2.MPI_SOURCE, 2,
                  this->cartTopology->GetComm(), &requestRecv2);
    }

    if (b1Sent) {
        MPI_Wait(&requestSend1, MPI_STATUS_IGNORE);
        MPI_Wait(&requestRecv1, MPI_STATUS_IGNORE);
    }
    if (b2Sent) {
        MPI_Wait(&requestSend2, MPI_STATUS_IGNORE);
        MPI_Wait(&requestRecv2, MPI_STATUS_IGNORE);
    }

    if (b1Sent) {
        this->b1 = this->b1Tmp;
        this->b1Tmp.clear();
    }

    if (b2Sent) {
        this->b2 = this->b2Tmp;
        this->b2Tmp.clear();
    }
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

int& P3BCA::getBufOwner(int i) {
    switch (i) {
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

std::tuple<uint64_t, uint64_t> P3BCA::SimulationStep() {
    this->b0 = this->simulation->GetDecomposition()->GetMyParticles();

    this->numShifts = 0;

#ifdef TESTS_3BMDA
    processed.clear();
#endif

#ifdef PROFILE_3BMDA
    this->hitrate = 0;
    uint64_t hitRateDivider = 0;
#endif

    uint64_t numBufferInteractions = 0;
    uint64_t numParticleInteractionsAcc = 0;

    // copy b0 to b1 and b2
    b1 = b0;
    b2 = b0;

    this->b1Owner = this->worldRank;
    this->b2Owner = this->worldRank;

    // 1D

    int nextSrcRankOuter1D;
    int nextDstRankOuter1D;

    int nextSrcRankInner1D;
    int nextDstRankInner1D;

    int diffOuter1D;
    int diffInner1D;

    int offsetVectorOuter1D = 0;
    int offsetVectorInner1D = 0;

    int myCoords1D;
    if (this->numDims == 1) {
        MPI_Cart_coords(this->cartTopology->GetComm(), this->worldRank, 1, &myCoords1D);
    }

    // 2D

    std::array<int, 2> nextSrcRankOuter2D;
    std::array<int, 2> nextDstRankOuter2D;

    std::array<int, 2> nextSrcRankInner2D;
    std::array<int, 2> nextDstRankInner2D;

    std::array<int, 2> diffOuter2D;
    std::array<int, 2> diffInner2D;

    std::array<int, 2> offsetVectorOuter2D = {0, 0};
    std::array<int, 2> offsetVectorInner2D = {0, 0};

    int myCoords2D[2];
    std::array<int, 2> myCoordsArray2D;
    if (this->numDims == 2) {
        MPI_Cart_coords(this->cartTopology->GetComm(), this->worldRank, 2, myCoords2D);
    }

    myCoordsArray2D[0] = myCoords2D[0];
    myCoordsArray2D[1] = myCoords2D[1];

    // 3D

    std::array<int, 3> nextSrcRankOuter3D;
    std::array<int, 3> nextDstRankOuter3D;

    std::array<int, 3> nextSrcRankInner3D;
    std::array<int, 3> nextDstRankInner3D;

    std::array<int, 3> diffOuter3D;
    std::array<int, 3> diffInner3D;

    std::array<int, 3> offsetVectorOuter3D = {0, 0, 0};
    std::array<int, 3> offsetVectorInner3D = {0, 0, 0};

    int myCoords3D[3];
    std::array<int, 3> myCoordsArray3D;
    if (this->numDims == 3) {
        MPI_Cart_coords(this->cartTopology->GetComm(), this->worldRank, 3, myCoords3D);
    }

    myCoordsArray3D[0] = myCoords3D[0];
    myCoordsArray3D[1] = myCoords3D[1];
    myCoordsArray3D[2] = myCoords3D[2];

    // edge case: #proc is 2
    for (int i2 = 0; i2 < (cartTopology->GetWorldSize() == 2 ? 1 : this->numSteps); i2++) {
        for (int i3 = i2; i3 < this->numSteps; i3++) {
#if defined(VLEVEL) && !defined(BENCHMARK_3BMDA) && !defined(TESTS_3BMDA) && VLEVEL > 0
            std::string message = "I'm proc " + std::to_string(simulation->GetTopology()->GetWorldRank()) +
                                  " and going to calculate interactions between buffers from proc (" +
                                  std::to_string(worldRank) + ", " + std::to_string(this->b1Owner) + ", " +
                                  std::to_string(this->b2Owner) + ")";
            MPIReporter::instance()->StoreMessage(this->simulation->GetTopology()->GetWorldRank(), message);
#endif

            std::tuple<uint64_t, uint64_t> numParticleInteractions =
                this->CalculateInteractions(this->b0, this->b1, this->b2, this->worldRank, this->b1Owner, this->b2Owner,
                                            this->cutoff, physicalDomainSize);
            numParticleInteractionsAcc += std::get<0>(numParticleInteractions);

            numBufferInteractions++;

#ifdef TESTS_3BMDA
            processed.push_back(Utility::Triplet(this->worldRank, this->b1Owner, this->b2Owner));
#endif

#ifdef PROFILE_3BMDA
            // only accumulate if there are possible particle interactions to avoid div by 0
            if (std::get<1>(numParticleInteractions) > 0) {
                this->hitrate += (double)std::get<0>(numParticleInteractions);
                hitRateDivider += std::get<1>(numParticleInteractions);
            }
#endif

            if (i3 < numSteps - 1) {
                if (this->numDims == 3) {
                    schedule3DHelper(i2, i3, myCoordsArray3D, nextSrcRankInner3D, nextDstRankInner3D, diffInner3D);
                } else if (this->numDims == 2) {
                    schedule2DHelper(i2, i3, myCoordsArray2D, nextSrcRankInner2D, nextDstRankInner2D, diffInner2D);
                } else {
                    schedule1DHelper(i3, myCoords1D, nextSrcRankInner1D, nextDstRankInner1D, diffInner1D);
                }

                if (i3 < this->numSteps - 1) {
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
                    if (this->numDims == 3) {
                        getBufOwner(2) = shiftLeft(this->b2, getBufOwner(2), nextSrcRankInner3D, nextDstRankInner3D,
                                                   offsetVectorInner3D, diffInner3D);
                    } else if (this->numDims == 2) {
                        getBufOwner(2) = shiftLeft(this->b2, getBufOwner(2), nextSrcRankInner2D, nextDstRankInner2D,
                                                   offsetVectorInner2D, diffInner2D);
                    } else {
                        getBufOwner(2) = shiftLeft(this->b2, getBufOwner(2), nextSrcRankInner1D, nextDstRankInner1D,
                                                   offsetVectorInner1D, diffInner1D);
                    }
                }
            }
        }

        // only shift if not the last step
        if (i2 < ((cartTopology->GetWorldSize() == 2 ? 1 : this->numSteps) - 1)) {
            if (this->numDims == 3) {
                schedule3D(i2 + 1, myCoordsArray3D, nextSrcRankOuter3D);
                calcDestFromSrc3D(myCoordsArray3D, nextSrcRankOuter3D, nextDstRankOuter3D);
                calcDiff3D(myCoordsArray3D, nextSrcRankOuter3D, diffOuter3D, i2 + 1);

                // also shift inner buffer
                schedule3D(i2 + 1, myCoordsArray3D, nextSrcRankInner3D);
                calcDestFromSrc3D(myCoordsArray3D, nextSrcRankInner3D, nextDstRankInner3D);
                calcDiff3D(myCoordsArray3D, nextSrcRankInner3D, diffInner3D, i2 + 1);

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

                getBufOwner(1) = shiftLeft(this->b1, getBufOwner(1), nextSrcRankOuter3D, nextDstRankOuter3D,
                                           offsetVectorOuter3D, diffOuter3D);

                getBufOwner(2) = shiftLeft(this->b2, getBufOwner(2), nextSrcRankInner3D, nextDstRankInner3D,
                                           offsetVectorInner3D, diffInner3D);

                offsetVectorInner3D = offsetVectorOuter3D;
            } else if (this->numDims == 2) {
                schedule2D(i2 + 1, myCoordsArray2D, nextSrcRankOuter2D);
                calcDestFromSrc2D(myCoordsArray2D, nextSrcRankOuter2D, nextDstRankOuter2D);
                calcDiff2D(myCoordsArray2D, nextSrcRankOuter2D, diffOuter2D, i2 + 1);

                // also shift inner buffer
                schedule2D(i2 + 1, myCoordsArray2D, nextSrcRankInner2D);
                calcDestFromSrc2D(myCoordsArray2D, nextSrcRankInner2D, nextDstRankInner2D);
                calcDiff2D(myCoordsArray2D, nextSrcRankInner2D, diffInner2D, i2 + 1);

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

                getBufOwner(1) = shiftLeft(this->b1, getBufOwner(1), nextSrcRankOuter2D, nextDstRankOuter2D,
                                           offsetVectorOuter2D, diffOuter2D);

                getBufOwner(2) = shiftLeft(this->b2, getBufOwner(2), nextSrcRankInner2D, nextDstRankInner2D,
                                           offsetVectorInner2D, diffInner2D);

                offsetVectorInner2D = offsetVectorOuter2D;
            } else {
                schedule1D(i2 + 1, myCoords1D, nextSrcRankOuter1D);
                calcDestFromSrc1D(myCoords1D, nextSrcRankOuter1D, nextDstRankOuter1D);
                calcDiff1D(myCoords1D, nextSrcRankOuter1D, diffOuter1D, i2 + 1);

                // also shift inner buffer
                schedule1D(i2 + 1, myCoords1D, nextSrcRankInner1D);
                calcDestFromSrc1D(myCoords1D, nextSrcRankInner1D, nextDstRankInner1D);
                calcDiff1D(myCoords1D, nextSrcRankInner1D, diffInner1D, i2 + 1);

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

                getBufOwner(1) = shiftLeft(this->b1, getBufOwner(1), nextSrcRankOuter1D, nextDstRankOuter1D,
                                           offsetVectorOuter1D, diffOuter1D);

                getBufOwner(2) = shiftLeft(this->b2, getBufOwner(2), nextSrcRankInner1D, nextDstRankInner1D,
                                           offsetVectorInner1D, diffInner1D);

                offsetVectorInner1D = offsetVectorOuter1D;
            }
        }
    }

    MPI_Barrier(this->simulation->GetTopology()->GetComm());

#ifdef PROFILE_3BMDA
    if (hitRateDivider > 0) {
        this->hitrate /= (double)hitRateDivider;
        this->hitrates.push_back(this->hitrate);
    }

#endif

    sendBackParticles();

    this->SumUpParticles(this->b0, this->b1, this->b2);

    this->simulation->GetDecomposition()->SetMyParticles(this->b0);

    return std::tuple(numBufferInteractions, numParticleInteractionsAcc);
}

std::array<int, 3> P3BCA::GetNumCutoffBoxes() { return std::array<int, 3>({this->nCbX, this->nCbY, this->nCbZ}); }

double P3BCA::GetCutoff() { return this->cutoff; }

std::tuple<uint64_t, uint64_t> P3BCA::SimulationStep([[maybe_unused]] ForceType forceType) {
    throw std::runtime_error("not implemented");
}