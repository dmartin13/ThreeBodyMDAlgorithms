#pragma once

#include "PairwisePotential.hpp"

class LennardJones final : public PairwisePotential {
private:
    const double epsilon = 1.0;
    const double sigma = 1.0;

public:
    LennardJones(double epsilon, double sigma);
    void CalculateForces(Utility::Particle &i, Utility::Particle &j) override;
};
