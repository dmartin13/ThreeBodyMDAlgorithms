#pragma once

#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <tuple>
#include <vector>

#include "DomainDecomposition.hpp"
#include "topology/CartTopology.hpp"

class RegularGridDecomposition final : public DomainDecomposition {
protected:
    CartRank cartRank;
    Eigen::Array3d localCellMin, localCellMax;
    std::shared_ptr<CartTopology> cartTopology;
    std::vector<Utility::Particle> sendToLeftNeighbor, sendToRightNeighbor;
    std::vector<Utility::Particle> recvFromLeftNeighbor, recvFromRightNeighbor;
    Eigen::Array3d localCellWidth;
    Eigen::Array3d physicalDomainSize;
    int dimX;
    int dimY;
    int dimZ;
    int numDims;

    void binParticles(std::vector<Utility::Particle>& particles);
    bool isInsideLocalCell(Utility::Particle& particle);
    std::tuple<Eigen::Array3d, Eigen::Array3d> getDomainMinMax(std::vector<Utility::Particle>& particles);
    void exchangeParticlesDim(int dim);
    void exchangeParticles();

public:
    RegularGridDecomposition();
    virtual ~RegularGridDecomposition();

    void Init(std::shared_ptr<Simulation> simulation) override;

    Eigen::Array3d GetCellSize();
    Eigen::Array3d GetPhysicalDomainSize();
    int GetDimX();
    int GetDimY();
    int GetDimZ();
};