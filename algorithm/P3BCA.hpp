#pragma once

#include "../decomposition/RegularGridDecomposition.hpp"
#include "../topology/CartTopology.hpp"
#include "Algorithm.hpp"

class P3BCA final : public Algorithm {
private:
    int worldRank;

    std::shared_ptr<CartTopology> cartTopology;
    std::vector<Utility::Particle>* b0;
    std::vector<Utility::Particle> b1;
    std::vector<Utility::Particle> b2;
    std::vector<Utility::Particle> tmpRecv;

    std::vector<Utility::Particle> b1Tmp;
    std::vector<Utility::Particle> b2Tmp;
    int b1Owner;
    int b2Owner;

    void calculateInteractions();
    int shift(std::vector<Utility::Particle>& buf, int dim, int dir, int owner);
    int shiftLeft(std::vector<Utility::Particle>& buf, int owner, std::array<int, 3>& nextSrcRank,
                  std::array<int, 3>& nextDstRank, std::array<int, 3>& offsetVector, std::array<int, 3>& diff);
    void shiftHelper(int i, std::array<int, 3>& myCartRank, std::array<int, 3>& src);
    void shiftHelper2(int i2, int& i3, std::array<int, 3>& cartRank, std::array<int, 3>& src,
                      std::array<int, 3>& dst, std::array<int, 3>& diff, int cutoffBorder);

    void sumUpParticles();
    void sendBackParticles();
    int& getBufOwner(int i);
    std::tuple<int, int> getInnerDirs(int i2);
    int periodicDistance(int x, int y, int dim);
    void moveForwardInLexManner();
    void calcDestFromSrc(std::array<int, 3>& myCartRank, std::array<int, 3>& src, std::array<int, 3>& dst);
    void calcDiff(std::array<int, 3>& cartRank, std::array<int, 3>& src, std::array<int, 3>& diff, int i);

    double cutoff;

    int numCutoffBoxes;

    int dim;

public:
    P3BCA(double cutoff);
    ~P3BCA();

    void Init(std::shared_ptr<Simulation> simulation) override;

    int SimulationStep() override;
    int SimulationStepNew();
    int SimulationStepNewNew();
    int GetNumCutoffBoxes();
};
