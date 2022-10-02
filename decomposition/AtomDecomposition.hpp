#pragma once

#include <iostream>
#include <vector>

#include "DomainDecomposition.hpp"

class AtomDecomposition final : public DomainDecomposition {
private:
    int worldRank;
    int worldSize;

public:
    AtomDecomposition();
    virtual ~AtomDecomposition();

    void Init(std::shared_ptr<Simulation> simulation) override;

    void Update(double dt, Eigen::Vector3d gForce) override;
    void UpdatePredictorStage(double dt) override;
};