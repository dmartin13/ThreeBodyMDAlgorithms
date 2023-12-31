#pragma once

#include <mpi.h>

#include <random>

#include "utility/structs.hpp"

namespace Thermostat {

    void addBrownianMotion(std::vector<Utility::Particle> &particles, const double targetTemperature);

    double calculateTemperature(std::vector<Utility::Particle> &particles);

    void apply(std::vector<Utility::Particle> &particles, const double targetTemperature,
               const double deltaTemperature);
}  // namespace Thermostat