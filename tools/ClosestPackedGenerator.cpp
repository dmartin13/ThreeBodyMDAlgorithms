#include "ClosestPackedGenerator.hpp"

ClosestPackedGenerator::ClosestPackedGenerator(int numParticles, const std::array<double, 3> &velocity,
                                               const std::array<double, 3> &boxLength,
                                               const std::array<double, 3> &bottomLeftCorner, double mass,
                                               uint_fast32_t seed0, uint_fast32_t seed1, double particleSpacing)
    : ParticleGenerator(numParticles, velocity, boxLength, bottomLeftCorner, mass, seed0, seed1),
      particleSpacing(particleSpacing)
{
    this->topRightCorner[0] = this->bottomLeftCorner[0] + this->boxLength[0];
    this->topRightCorner[1] = this->bottomLeftCorner[1] + this->boxLength[1];
    this->topRightCorner[2] = this->bottomLeftCorner[2] + this->boxLength[2];
}

ClosestPackedGenerator::~ClosestPackedGenerator() {}

void ClosestPackedGenerator::Generate()
{
    const double spacingRow = particleSpacing * sqrt(3. / 4.);
    const double spacingLayer = particleSpacing * sqrt(2. / 3.);
    const double xOffset = particleSpacing * 1. / 2.;
    const double yOffset = particleSpacing * sqrt(1. / 12.);

    bool evenLayer = true;
    bool evenRow = true;

    int id = 0;
    for (double z = bottomLeftCorner[2]; z < topRightCorner[2]; z += spacingLayer) {
        double starty = evenLayer ? bottomLeftCorner[1] : bottomLeftCorner[1] + yOffset;
        for (double y = starty; y < topRightCorner[1]; y += spacingRow) {
            double startx = evenRow ? bottomLeftCorner[0] : bottomLeftCorner[0] + xOffset;
            for (double x = startx; x < topRightCorner[0]; x += particleSpacing) {
                std::tuple<int, double, double, double, double, double, double, double, double, double, double>
                    positions = std::make_tuple(id++, x, y, z, velocity[0], velocity[1], velocity[2], 0, 0, 0, mass);

                this->particles.push_back(positions);
            }
            evenRow = !evenRow;
        }
        evenLayer = !evenLayer;
    }
}
