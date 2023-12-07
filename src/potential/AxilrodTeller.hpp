#pragma once

#include "TriwisePotential.hpp"

class AxilrodTeller final : public TriwisePotential {
private:
    const double nu = 1.0;

public:
    AxilrodTeller(double nu);
    void CalculateForces(Utility::Particle &i, Utility::Particle &j, Utility::Particle &k) override;
};
