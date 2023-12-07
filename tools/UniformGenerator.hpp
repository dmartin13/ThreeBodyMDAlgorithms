#pragma once

#include "ParticleGenerator.hpp"

class UniformGenerator : public ParticleGenerator {
private:
    /* data */
public:
    UniformGenerator(int numParticles, const std::array<double, 3> &velocity, const std::array<double, 3> &boxLength,
                     const std::array<double, 3> &bottomLeftCorner, double mass, uint_fast32_t seed0,
                     uint_fast32_t seed1);
    ~UniformGenerator();

    void Generate() override;
};
