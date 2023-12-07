#pragma once

#include "ParticleGenerator.hpp"

class ClusteredGaussGenerator : public ParticleGenerator {
private:
    const std::array<double, 3> &distributionMean;
    const std::array<double, 3> &distributionStdDev;
    int numClusters;
public:
    ClusteredGaussGenerator(int numParticles, const std::array<double, 3> &velocity,
                            const std::array<double, 3> &boxLength, const std::array<double, 3> &bottomLeftCorner,
                            double mass, uint_fast32_t seed0, uint_fast32_t seed1,
                            const std::array<double, 3> &distributionMean,
                            const std::array<double, 3> &distributionStdDev, int numClusters);
    ~ClusteredGaussGenerator();

    void Generate() override;
};
