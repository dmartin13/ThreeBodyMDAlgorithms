#include "AxilrodTeller.hpp"

AxilrodTeller::AxilrodTeller(double nu) : nu(nu) {}

void AxilrodTeller::CalculateForces(Utility::Particle &i, Utility::Particle &j, Utility::Particle &k) {
    // we assume forces are set to 0 for each particle

    // sides of triangle
    Eigen::Vector3d a = (j.GetR() - i.GetR());
    Eigen::Vector3d b = (i.GetR() - k.GetR());
    Eigen::Vector3d c = (j.GetR() - k.GetR());

    // Explanation of notation: d = distance, number = power, a = (i,j), b = (i, k), c = (j, k)

    // distances and powers of distances
    double d1a = a.norm();
    double d1b = b.norm();
    double d1c = c.norm();

    double d2a = d1a * d1a;
    double d2b = d1b * d1b;
    double d2c = d1c * d1c;

    double d3a = d2a * d1a;
    double d3b = d2b * d1b;
    double d3c = d2c * d1c;

    double d4a = d3a * d1a;
    double d4b = d3b * d1b;
    double d4c = d3c * d1c;

    double d5a = d4a * d1a;
    double d5b = d4b * d1b;
    double d5c = d4c * d1c;

    double d6a = d5a * d1a;
    double d6b = d5b * d1b;
    double d6c = d5c * d1c;

    double dXa = i.pX - j.pX;
    double dYa = i.pY - j.pY;
    double dZa = i.pZ - j.pZ;

    double dXb = i.pX - k.pX;
    double dYb = i.pY - k.pY;
    double dZb = i.pZ - k.pZ;

    double dXc = j.pX - k.pX;
    double dYc = j.pY - k.pY;
    double dZc = j.pZ - k.pZ;

    double dVdRa, dVdRb, dVdRc;

    // this is the gradient of the axilrodteller potential
    dVdRa = (3. * nu / (8. * d1a)) *
            (-8. / (d4a * d3b * d3c) - 1. / (d5b * d5c) + 5. * d1b / (d6a * d5c) + 5. * d1c / (d6a * d5b) -
             1. / (d2a * d3b * d5c) - 1. / (d2a * d5b * d3c) - 3. / (d4a * d1b * d5c) - 3. / (d4a * d5b * d1c) -
             5. / (d6a * d1b * d3c) - 5. / (d6a * d3b * d1c) + 6. / (d4a * d3b * d3c));

    dVdRb = (3. * nu / (8. * d1b)) *
            (-8. / (d4b * d3a * d3c) - 1. / (d5a * d5c) + 5. * d1a / (d6b * d5c) + 5. * d1c / (d6b * d5a) -
             1. / (d2b * d3a * d5c) - 1. / (d2b * d5a * d3c) - 3. / (d4b * d1a * d5c) - 3. / (d4b * d5a * d1c) -
             5. / (d6b * d1a * d3c) - 5. / (d6b * d3a * d1c) + 6. / (d4b * d3a * d3c));

    dVdRc = (3. * nu / (8. * d1c)) *
            (-8. / (d4c * d3b * d3a) - 1. / (d5b * d5a) + 5. * d1b / (d6c * d5a) + 5. * d1a / (d6c * d5b) -
             1. / (d2c * d3b * d5a) - 1. / (d2c * d5b * d3a) - 3. / (d4c * d1b * d5a) - 3. / (d4c * d5b * d1a) -
             5. / (d6c * d1b * d3a) - 5. / (d6c * d3b * d1a) + 6. / (d4c * d3b * d3a));

    auto forceIX = -dXa * dVdRa - dXb * dVdRb;
    auto forceIY = -dYa * dVdRa - dYb * dVdRb;
    auto forceIZ = -dZa * dVdRa - dZb * dVdRb;

    auto forceJX = -dXa * (-dVdRa) - dXc * dVdRc;
    auto forceJY = -dYa * (-dVdRa) - dYc * dVdRc;
    auto forceJZ = -dZa * (-dVdRa) - dZc * dVdRc;

    auto forceKX = -dXb * (-dVdRb) - dXc * (-dVdRc);
    auto forceKY = -dYb * (-dVdRb) - dYc * (-dVdRc);
    auto forceKZ = -dZb * (-dVdRb) - dZc * (-dVdRc);

    i.f1X += forceIX;
    i.f1Y += forceIY;
    i.f1Z += forceIZ;

    j.f1X += forceJX;
    j.f1Y += forceJY;
    j.f1Z += forceJZ;

    k.f1X += forceKX;
    k.f1Y += forceKY;
    k.f1Z += forceKZ;

    // potential energy
    const auto upot = nu * ((1.0 / (d3a * d3b * d3b)) +
                            (3 * ((-d2a + d2b + d2c) * (d2a - d2b + d2c) * (d2a + d2b - d2c))) / (8 * d5a * d5b * d5c));
    potentialEnergy += upot;

    // virial
    virial +=
        (Eigen::Array3d{forceIX, forceIY, forceIZ} * i.GetR() + Eigen::Array3d{forceJX, forceJY, forceJZ} * j.GetR() +
         Eigen::Array3d{forceKX, forceKY, forceKZ} * k.GetR());
}