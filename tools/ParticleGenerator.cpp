#include "ParticleGenerator.hpp"

ParticleGenerator::ParticleGenerator(int numParticles, const std::array<double, 3> &velocity,
                                     const std::array<double, 3> &boxLength,
                                     const std::array<double, 3> &bottomLeftCorner, double mass, uint_fast32_t seed0,
                                     uint_fast32_t seed1)
    : numParticles(numParticles), velocity(velocity), boxLength(boxLength), bottomLeftCorner(bottomLeftCorner),
      mass(mass), seed0(seed0), seed1(seed1)
{}

ParticleGenerator::~ParticleGenerator() {}

std::vector<std::tuple<int, double, double, double, double, double, double, double, double, double, double>>
ParticleGenerator::GetParticles()
{
    return this->particles;
}

ParticleGenerator::Generator ParticleGenerator::Str2Gen(std::string str)
{
    if (str.compare("closestpacked") == 0) {
        return Generator::ClosestPacked;
    } else if (str.compare("clusteredgauss") == 0) {
        return Generator::ClusteredGauss;
    } else if (str.compare("gauss") == 0) {
        return Generator::Gauss;
    } else if (str.compare("grid") == 0) {
        return Generator::Grid;
    } else if (str.compare("uniform") == 0) {
        return Generator::Uniform;
    } else {
        return Generator::Uniform;
    }
}