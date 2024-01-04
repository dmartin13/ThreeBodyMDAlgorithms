#pragma once

#include "ParticleGenerator.hpp"

class GaussGenerator : public ParticleGenerator {
private:
    const std::array<double, 3> &distributionMean;
    const std::array<double, 3> &distributionStdDev;

public:
    GaussGenerator(int numParticles, const std::array<double, 3> &velocity, const std::array<double, 3> &boxLength,
                   const std::array<double, 3> &bottomLeftCorner, double mass, uint_fast32_t seed0, uint_fast32_t seed1,
                   const std::array<double, 3> &distributionMean, const std::array<double, 3> &distributionStdDev);
    ~GaussGenerator();

    void Generate() override;
};