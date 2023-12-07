#pragma once

#include "PairwisePotential.hpp"

class LennardJones final : public PairwisePotential {
private:
    const double _epsilon = 1.0;
    const double _sigma = 1.0;
    const double _sigmaSquared = 1.0;
    const double _epsilon24 = 24.0;
    const double _cutoffSquared = 0;

public:
    LennardJones(double epsilon, double sigma);
    void CalculateForces(Utility::Particle &i, Utility::Particle &j) override;
};
