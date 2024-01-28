#include "TimeIntegration.hpp"

namespace TimeIntegration {
    // this based on the integration schema in Griebel
    void calculatePositionsAndResetForces(std::vector<Utility::Particle>& particles, const double deltaT,
                                          const std::array<double, 3> globalForce, ForceType forceTypeToUse) {
        for (auto& p : particles) {
            if (p.isDummy) {
                continue;
            }

            // read forces
            auto fX = p.f0X;
            auto fY = p.f0Y;
            auto fZ = p.f0Z;
            if (forceTypeToUse == ForceType::TwoAndThreeBody) {
                fX += p.f1X;
                fY += p.f1Y;
                fZ += p.f1Z;
            }

            // set old force
            p.oldF0X = p.f0X;
            p.oldF0Y = p.f0Y;
            p.oldF0Z = p.f0Z;
            if (forceTypeToUse == ForceType::TwoAndThreeBody) {
                p.oldF0X += p.f1X;
                p.oldF0Y += p.f1Y;
                p.oldF0Z += p.f1Z;
            }

            // reset force to the global Force
            p.f0X = globalForce[0];
            p.f0Y = globalForce[1];
            p.f0Z = globalForce[2];
            // if (forceTypeToUse == ForceType::TwoAndThreeBody) {
            p.f1X = 0;
            p.f1Y = 0;
            p.f1Z = 0;
            //}

            // velocity
            auto vX = p.vX * deltaT;
            auto vY = p.vY * deltaT;
            auto vZ = p.vZ * deltaT;

            fX *= (deltaT * deltaT / (2 * p.mass));
            fY *= (deltaT * deltaT / (2 * p.mass));
            fZ *= (deltaT * deltaT / (2 * p.mass));

            // calcuate displacement
            const auto displacementX = vX + fX;
            const auto displacementY = vY + fY;
            const auto displacementZ = vZ + fZ;

            // update positions
            p.pX += displacementX;
            p.pY += displacementY;
            p.pZ += displacementZ;
        }
    }

    void calculateVelocities(std::vector<Utility::Particle>& particles, const double deltaT, const bool outerRespaStep,
                             const int respaStepSize) {
        for (auto& p : particles) {
            if (p.isDummy) {
                continue;
            }
            if (respaStepSize == -1) {
                // respa is not active -> use both two and threebody force to update the velocities
                const auto fX = p.f0X + p.f1X;
                const auto fY = p.f0Y + p.f1Y;
                const auto fZ = p.f0Z + p.f1Z;

                // Is this correct or do we need p.oldF0X + p.oldF1X here? -> Is correct since for non respa simulations
                // we write both p.f0X + p.f1X into p.oldF0X. See calculatePositionsAndResetForces above
                const auto oldForceX = p.oldF0X;
                const auto oldForceY = p.oldF0Y;
                const auto oldForceZ = p.oldF0Z;

                const auto changeInVelX = (fX + oldForceX) * (deltaT / (2 * p.mass));
                const auto changeInVelY = (fY + oldForceY) * (deltaT / (2 * p.mass));
                const auto changeInVelZ = (fZ + oldForceZ) * (deltaT / (2 * p.mass));

                p.vX += changeInVelX;
                p.vY += changeInVelY;
                p.vZ += changeInVelZ;
            } else {
                if (not outerRespaStep) {
                    // inner respa step where we use the two-body force to update velocities
                    const auto fX = p.f0X;
                    const auto fY = p.f0Y;
                    const auto fZ = p.f0Z;

                    const auto oldForceX = p.oldF0X;
                    const auto oldForceY = p.oldF0Y;
                    const auto oldForceZ = p.oldF0Z;

                    const auto changeInVelX = (fX + oldForceX) * (deltaT / (2 * p.mass));
                    const auto changeInVelY = (fY + oldForceY) * (deltaT / (2 * p.mass));
                    const auto changeInVelZ = (fZ + oldForceZ) * (deltaT / (2 * p.mass));

                    p.vX += changeInVelX;
                    p.vY += changeInVelY;
                    p.vZ += changeInVelZ;
                } else {
                    // outer respa step where we use three-body force to update velocity
                    const auto fX = p.f1X;
                    const auto fY = p.f1Y;
                    const auto fZ = p.f1Z;

                    const auto changeInVelX = (fX * (1.0 / p.mass)) * 0.5 * deltaT * static_cast<double>(respaStepSize);
                    const auto changeInVelY = (fY * (1.0 / p.mass)) * 0.5 * deltaT * static_cast<double>(respaStepSize);
                    const auto changeInVelZ = (fZ * (1.0 / p.mass)) * 0.5 * deltaT * static_cast<double>(respaStepSize);

                    p.vX += changeInVelX;
                    p.vY += changeInVelY;
                    p.vZ += changeInVelZ;
                }
            }
        }
    }
}  // namespace TimeIntegration