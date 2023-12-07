#include "GaussGenerator.hpp"

GaussGenerator::GaussGenerator(int numParticles, const std::array<double, 3> &velocity,
                               const std::array<double, 3> &boxLength, const std::array<double, 3> &bottomLeftCorner,
                               double mass, uint_fast32_t seed0, uint_fast32_t seed1,
                               const std::array<double, 3> &distributionMean,
                               const std::array<double, 3> &distributionStdDev)
    : ParticleGenerator(numParticles, velocity, boxLength, bottomLeftCorner, mass, seed0, seed1),
      distributionMean(distributionMean), distributionStdDev(distributionStdDev)
{}

GaussGenerator::~GaussGenerator() {}

void GaussGenerator::Generate()
{
    std::default_random_engine generator(this->seed0);
    std::array<std::normal_distribution<double>, 3> distributions = {
        std::normal_distribution<double>{distributionMean[0], distributionStdDev[0]},
        std::normal_distribution<double>{distributionMean[1], distributionStdDev[1]},
        std::normal_distribution<double>{distributionMean[2], distributionStdDev[2]}};

    for (int i = 0; i < numParticles; ++i) {
        std::tuple<int, double, double, double, double, double, double, double, double, double, double> positions =
            std::make_tuple(i, bottomLeftCorner[0] + distributions[0](generator),
                            bottomLeftCorner[1] + distributions[1](generator),
                            bottomLeftCorner[2] + distributions[2](generator), velocity[0], velocity[1], velocity[2], 0,
                            0, 0, mass);
        this->particles.push_back(positions);
    }
}