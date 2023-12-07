#include "ClusteredGaussGenerator.hpp"

ClusteredGaussGenerator::ClusteredGaussGenerator(int numParticles, const std::array<double, 3> &velocity,
                                                 const std::array<double, 3> &boxLength,
                                                 const std::array<double, 3> &bottomLeftCorner, double mass,
                                                 uint_fast32_t seed0, uint_fast32_t seed1,
                                                 const std::array<double, 3> &distributionMean,
                                                 const std::array<double, 3> &distributionStdDev, int numClusters)
    : ParticleGenerator(numParticles, velocity, boxLength, bottomLeftCorner, mass, seed0, seed1),
      distributionMean(distributionMean), distributionStdDev(distributionStdDev), numClusters(numClusters)
{}

ClusteredGaussGenerator::~ClusteredGaussGenerator() {}

void ClusteredGaussGenerator::Generate()
{
    // create random points in box
    std::vector<std::tuple<double, double, double>> points;

    std::mt19937 randomNumberEngine(seed1);
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    for (int i = 0; i < numClusters; ++i) {
        std::tuple<double, double, double> positions =

            std::make_tuple(bottomLeftCorner[0] + distribution(randomNumberEngine) * boxLength[0],
                            bottomLeftCorner[1] + distribution(randomNumberEngine) * boxLength[1],
                            bottomLeftCorner[2] + distribution(randomNumberEngine) * boxLength[2]);
        points.push_back(positions);
    }

    std::default_random_engine generator(this->seed0);
    std::array<std::normal_distribution<double>, 3> distributions = {
        std::normal_distribution<double>{distributionMean[0], distributionStdDev[0]},
        std::normal_distribution<double>{distributionMean[1], distributionStdDev[1]},
        std::normal_distribution<double>{distributionMean[2], distributionStdDev[2]}};

    int id = 0;
    for (std::tuple<double, double, double> &p : points) {
        for (int i = 0; i < numParticles / numClusters; ++i) {
            std::tuple<int, double, double, double, double, double, double, double, double, double, double> positions =
                std::make_tuple(
                    id++, std::get<0>(p) + distributions[0](generator), std::get<1>(p) + distributions[1](generator),
                    std::get<2>(p) + distributions[2](generator), velocity[0], velocity[1], velocity[2], 0, 0, 0, mass);
            this->particles.push_back(positions);
        }
    }
}