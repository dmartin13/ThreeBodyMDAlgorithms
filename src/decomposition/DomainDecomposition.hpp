
#pragma once

#include <memory>
#include <vector>

#include "fwd.hpp"
#include "simulation/Simulation.hpp"
#include "utility/structs.hpp"

class DomainDecomposition {
protected:
    std::shared_ptr<Simulation> simulation;
    std::vector<Utility::Particle> myParticles;

    void updateMyParticles(double dt, Eigen::Vector3d gForce);
    void updateMyParticlesPredictorStage(double dt);

public:
    DomainDecomposition();
    virtual ~DomainDecomposition();

    virtual void Init(std::shared_ptr<Simulation> simulation);

    virtual void Update(double dt, Eigen::Vector3d gForce) = 0;
    virtual void UpdatePredictorStage(double dt) = 0;
    void ResetForces();
    std::vector<Utility::Particle> GetMyParticles();
    void SetMyParticles(std::vector<Utility::Particle> &particles);
    int GetNumOfMyParticles();
};