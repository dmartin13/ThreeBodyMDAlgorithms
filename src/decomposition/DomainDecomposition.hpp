
#pragma once

#include <memory>
#include <vector>

#include "fwd.hpp"
#include "simulation/Simulation.hpp"
#include "simulation/timeintegration/TimeIntegration.hpp"
#include "utility/structs.hpp"

class DomainDecomposition {
protected:
    std::shared_ptr<Simulation> simulation;
    std::vector<Utility::Particle> myParticles;

public:
    DomainDecomposition();
    virtual ~DomainDecomposition();

    virtual void Init(std::shared_ptr<Simulation> simulation);

    void UpdatePositions(double dt, Eigen::Vector3d gForce, ForceType forceTypeToUse);
    void UpdateVelocities(double dt, ForceType forceType, size_t respaStepSize);
    std::vector<Utility::Particle>& GetMyParticles();
    void SetMyParticles(std::vector<Utility::Particle>& particles);
    int GetNumOfMyParticles();
};