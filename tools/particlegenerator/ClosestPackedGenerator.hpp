#pragma once

#include "ParticleGenerator.hpp"

class ClosestPackedGenerator : public ParticleGenerator {
private:
    double particleSpacing;
    std::array<double, 3> topRightCorner;

public:
    ClosestPackedGenerator(int numParticles, const std::array<double, 3> &velocity,
                           const std::array<double, 3> &boxLength, const std::array<double, 3> &bottomLeftCorner,
                           double mass, uint_fast32_t seed0, uint_fast32_t seed1, double particleSpacing);
    ~ClosestPackedGenerator();

    void Generate() override;
};