#pragma once

#include "ParticleGenerator.hpp"

class GridGenerator : public ParticleGenerator {
private:
    std::array<size_t, 3> particlesPerDim;
    double particleSpacing;

public:
    GridGenerator(int numParticles, const std::array<double, 3> &velocity, const std::array<double, 3> &boxLength,
                  const std::array<double, 3> &bottomLeftCorner, double mass, uint_fast32_t seed0, uint_fast32_t seed1,
                  std::array<size_t, 3> &particlesPerDim, double particleSpacing);
    ~GridGenerator();

    void Generate() override;
};