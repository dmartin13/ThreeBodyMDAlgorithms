#pragma once

#include <array>

#include "utility/enums.hpp"
#include "utility/structs.hpp"

namespace TimeIntegration {
    // this based on the integration schema in Griebel
    void calculatePositionsAndResetForces(std::vector<Utility::Particle>& particles, const double deltaT,
                                          const std::array<double, 3> globalForce, ForceType forceTypeToUse);

    void calculateVelocities(std::vector<Utility::Particle>& particles, const double deltaT, const bool outerRespaStep,
                             const int respaStepSize);
}  // namespace TimeIntegration