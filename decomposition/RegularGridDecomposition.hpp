#pragma once

#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <tuple>
#include <vector>

#include "../topology/CartTopology.hpp"
#include "DomainDecomposition.hpp"

class RegularGridDecomposition final : public DomainDecomposition {
protected:
    CartRank cartRank;
    Eigen::Array3d localCellMin, localCellMax;
    std::shared_ptr<CartTopology> cartTopology;
    std::vector<Utility::Particle> sendToLeftNeighbor, sendToRightNeighbor;
    std::vector<Utility::Particle> recvFromLeftNeighbor, recvFromRightNeighbor;
    Eigen::Array3d localCellWidth;
    // double physicalDomainSize;
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

    void Update(double dt, Eigen::Vector3d gForce) override;

    Eigen::Array3d GetCellSize();
    int GetDimX();
    int GetDimY();
    int GetDimZ();
};