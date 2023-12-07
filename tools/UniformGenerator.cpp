#include "UniformGenerator.hpp"

UniformGenerator::UniformGenerator(int numParticles, const std::array<double, 3> &velocity,
                                   const std::array<double, 3> &boxLength,
                                   const std::array<double, 3> &bottomLeftCorner, double mass, uint_fast32_t seed0,
                                   uint_fast32_t seed1)
    : ParticleGenerator(numParticles, velocity, boxLength, bottomLeftCorner, mass, seed0, seed1)
{}

UniformGenerator::~UniformGenerator() {}

void UniformGenerator::Generate()
{
    // Set up random number generation
    // std::random_device randomDevice;
    std::mt19937 randomNumberEngine(seed0);
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    for (int i = 0; i < numParticles; ++i) {
        std::tuple<int, double, double, double, double, double, double, double, double, double, double> positions =

            std::make_tuple(i, bottomLeftCorner[0] + distribution(randomNumberEngine) * boxLength[0],
                            bottomLeftCorner[1] + distribution(randomNumberEngine) * boxLength[1],
                            bottomLeftCorner[2] + distribution(randomNumberEngine) * boxLength[2], velocity[0],
                            velocity[1], velocity[2], 0, 0, 0, mass);
        this->particles.push_back(positions);
    }
}