#pragma once

#include "decomposition/RegularGridDecomposition.hpp"
#include "topology/CartTopology.hpp"
#include "Algorithm.hpp"

class P3BCA final : public Algorithm {
private:
    double cutoff;
    int worldRank;
    int numDims;
    int dimX;
    int dimY;
    int dimZ;
    int nCbX;
    int nCbY;
    int nCbZ;
    int numSteps;
    Eigen::Array3d physicalDomainSize;

    int b1Owner;
    int b2Owner;

    std::shared_ptr<CartTopology> cartTopology;
    std::vector<Utility::Particle> b0;
    std::vector<Utility::Particle> b1;
    std::vector<Utility::Particle> b2;

    std::vector<Utility::Particle> tmpRecv;
    std::vector<Utility::Particle> b1Tmp;
    std::vector<Utility::Particle> b2Tmp;

    void calcSteps(int dimension);

    int shiftLeft(std::vector<Utility::Particle>& buf, int owner, int& nextSrcRank, int& nextDstRank, int&
       offsetVector, int& diff);
    int shiftLeft(std::vector<Utility::Particle>& buf, int owner, std::array<int, 2>& nextSrcRank,
                  std::array<int, 2>& nextDstRank, std::array<int, 2>& offsetVector, std::array<int, 2>& diff);
    int shiftLeft(std::vector<Utility::Particle>& buf, int owner, std::array<int, 3>& nextSrcRank,
                  std::array<int, 3>& nextDstRank, std::array<int, 3>& offsetVector, std::array<int, 3>& diff);
    int mpiShift(std::vector<Utility::Particle>& buf, int owner, int src, int dst);

    void sendBackParticles();

    int& getBufOwner(int i);

    // 1D
    
    void schedule1D(int i, int& myCartRank, int& src);
    void calcDestFromSrc1D(int& myCartRank, int& src, int& dst);
    void schedule1DHelper(int i3, int& cartRank, int& src, int& dst, int& diff);
    void calcDiff1D(int& cartRank, int& src, int& diff, int i);
    void handleOffsetVector1D(int& nextSrcRank, int& nextDstRank, int& offsetVector, int& diff,
                              int& coordsSrc, int& coordsDst);
    

    // 2D
    void schedule2D(int i, std::array<int, 2>& myCartRank, std::array<int, 2>& src);
    void calcDestFromSrc2D(std::array<int, 2>& myCartRank, std::array<int, 2>& src, std::array<int, 2>& dst);
    void schedule2DHelper(int i2, int& i3, std::array<int, 2>& cartRank, std::array<int, 2>& src,
                          std::array<int, 2>& dst, std::array<int, 2>& diff);
    void calcDiff2D(std::array<int, 2>& cartRank, std::array<int, 2>& src, std::array<int, 2>& diff, int i);
    void handleOffsetVector2D(std::array<int, 2>& nextSrcRank, std::array<int, 2>& nextDstRank,
                              std::array<int, 2>& offsetVector, std::array<int, 2>& diff, std::array<int, 2>& coordsSrc,
                              std::array<int, 2>& coordsDst);

    // 3D
    void schedule3D(int i, std::array<int, 3>& myCartRank, std::array<int, 3>& src);
    void calcDestFromSrc3D(std::array<int, 3>& myCartRank, std::array<int, 3>& src, std::array<int, 3>& dst);
    void schedule3DHelper(int i2, int& i3, std::array<int, 3>& cartRank, std::array<int, 3>& src,
                          std::array<int, 3>& dst, std::array<int, 3>& diff);
    void calcDiff3D(std::array<int, 3>& cartRank, std::array<int, 3>& src, std::array<int, 3>& diff, int i);
    void handleOffsetVector3D(std::array<int, 3>& nextSrcRank, std::array<int, 3>& nextDstRank,
                              std::array<int, 3>& offsetVector, std::array<int, 3>& diff, std::array<int, 3>& coordsSrc,
                              std::array<int, 3>& coordsDst);

public:
    P3BCA(double cutoff);
    virtual ~P3BCA();

    void Init(std::shared_ptr<Simulation> simulation) override;

    std::tuple<uint64_t, uint64_t> SimulationStep() override;
    int SimulationStepNew();
    int SimulationStepNewNew();
    std::array<int, 3> GetNumCutoffBoxes();
    double GetCutoff();
};