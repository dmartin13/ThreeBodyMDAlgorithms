#include "Thermostat.hpp"

namespace Thermostat {

    void addBrownianMotion(std::vector<Utility::Particle> &particles, const double targetTemperature) {
        // this assumes that all particles have the same mass
        const auto translationalVelocityScale = std::sqrt(targetTemperature / particles[0].mass);

        // we use a constant seed for repeatability.
        std::default_random_engine randomEngine(42);
        std::normal_distribution<double> normalDistribution{0, 1};

        for (auto &p : particles) {
            if (p.isDummy) {
                continue;
            }
            p.vX += normalDistribution(randomEngine) * translationalVelocityScale;
            p.vY += normalDistribution(randomEngine) * translationalVelocityScale;
            p.vZ += normalDistribution(randomEngine) * translationalVelocityScale;
        }
    }

    double calculateTemperature(std::vector<Utility::Particle> &particles) {
        double result = 0;

        for (const auto &p : particles) {
            if (p.isDummy) {
                continue;
            }
            // kinetic energy * 2
            result += p.mass * (p.vX * p.vX + p.vY * p.vY + p.vZ * p.vZ);
        }

        // get number of particles
        size_t numParticles =
            std::count_if(particles.begin(), particles.end(), [](const auto &p) { return (not p.isDummy); });

        MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &numParticles, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);

        // divide by number of particles * 3 DoF
        result /= static_cast<double>(numParticles) * 3.0;

        return result;
    }

    void apply(std::vector<Utility::Particle> &particles, const double targetTemperature,
               const double deltaTemperature) {
        const auto currentTemperature = calculateTemperature(particles);

        // make sure we work with a positive delta
        const double absoluteDeltaTemperature = std::abs(deltaTemperature);

        // calculate scaling factor
        const auto immediateTargetTemperature =
            currentTemperature < targetTemperature
                ? std::min(currentTemperature + absoluteDeltaTemperature, targetTemperature)
                : std::max(currentTemperature - absoluteDeltaTemperature, targetTemperature);

        const auto scalingFactor = std::sqrt(immediateTargetTemperature / currentTemperature);

        for (auto &p : particles) {
            if (p.isDummy) {
                continue;
            }
            p.vX *= scalingFactor;
            p.vY *= scalingFactor;
            p.vZ *= scalingFactor;
        }
    }
}  // namespace Thermostat