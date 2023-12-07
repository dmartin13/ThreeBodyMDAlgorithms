#include "GridGenerator.hpp"

GridGenerator::GridGenerator(int numParticles, const std::array<double, 3> &velocity,
                             const std::array<double, 3> &boxLength, const std::array<double, 3> &bottomLeftCorner,
                             double mass, uint_fast32_t seed0, uint_fast32_t seed1,
                             std::array<size_t, 3> &particlesPerDim, double particleSpacing)
    : ParticleGenerator(numParticles, velocity, boxLength, bottomLeftCorner, mass, seed0, seed1),
      particlesPerDim(particlesPerDim), particleSpacing(particleSpacing)
{}

GridGenerator::~GridGenerator() {}

void GridGenerator::Generate()
{
    int id = 0;
    for (unsigned long z = 0; z < particlesPerDim[2]; ++z) {
        for (unsigned long y = 0; y < particlesPerDim[1]; ++y) {
            for (unsigned long x = 0; x < particlesPerDim[0]; ++x) {
                std::tuple<int, double, double, double, double, double, double, double, double, double, double>
                    positions =

                        std::make_tuple(id++, bottomLeftCorner[0] + static_cast<double>(x) * particleSpacing,
                                        bottomLeftCorner[1] + static_cast<double>(y) * particleSpacing,
                                        bottomLeftCorner[2] + static_cast<double>(z) * particleSpacing, velocity[0],
                                        velocity[1], velocity[2], 0, 0, 0, mass);
                this->particles.push_back(positions);
            }
        }
    }
}
