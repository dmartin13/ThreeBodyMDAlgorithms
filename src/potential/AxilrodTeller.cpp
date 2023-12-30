#include "AxilrodTeller.hpp"

AxilrodTeller::AxilrodTeller(double nu) : nu(nu) {}

void AxilrodTeller::CalculateForces(Utility::Particle &i, Utility::Particle &j, Utility::Particle &k) {
    // we assume forces are set to 0 for each particle

    // sides of triangle
    Eigen::Vector3d a = (j.GetR() - i.GetR());
    Eigen::Vector3d b = (i.GetR() - k.GetR());
    Eigen::Vector3d c = (j.GetR() - k.GetR());

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
    dVdRa = (3. * this->nu / (8. * d1a)) *
            (-8. / (d4a * d3b * d3c) - 1. / (d5b * d5c) + 5. * d1b / (d6a * d5c) + 5. * d1c / (d6a * d5b) -
             1. / (d2a * d3b * d5c) - 1. / (d2a * d5b * d3c) - 3. / (d4a * d1b * d5c) - 3. / (d4a * d5b * d1c) -
             5. / (d6a * d1b * d3c) - 5. / (d6a * d3b * d1c) + 6. / (d4a * d3b * d3c));

    dVdRb = (3. * this->nu / (8. * d1b)) *
            (-8. / (d4b * d3a * d3c) - 1. / (d5a * d5c) + 5. * d1a / (d6b * d5c) + 5. * d1c / (d6b * d5a) -
             1. / (d2b * d3a * d5c) - 1. / (d2b * d5a * d3c) - 3. / (d4b * d1a * d5c) - 3. / (d4b * d5a * d1c) -
             5. / (d6b * d1a * d3c) - 5. / (d6b * d3a * d1c) + 6. / (d4b * d3a * d3c));

    dVdRc = (3. * this->nu / (8. * d1c)) *
            (-8. / (d4c * d3b * d3a) - 1. / (d5b * d5a) + 5. * d1b / (d6c * d5a) + 5. * d1a / (d6c * d5b) -
             1. / (d2c * d3b * d5a) - 1. / (d2c * d5b * d3a) - 3. / (d4c * d1b * d5a) - 3. / (d4c * d5b * d1a) -
             5. / (d6c * d1b * d3a) - 5. / (d6c * d3b * d1a) + 6. / (d4c * d3b * d3a));

    i.f1X = i.f1X - dXa * dVdRa - dXb * dVdRb;
    i.f1Y = i.f1Y - dYa * dVdRa - dYb * dVdRb;
    i.f1Z = i.f1Z - dZa * dVdRa - dZb * dVdRb;

    j.f1X = j.f1X - dXa * (-dVdRa) - dXc * dVdRc;
    j.f1Y = j.f1Y - dYa * (-dVdRa) - dYc * dVdRc;
    j.f1Z = j.f1Z - dZa * (-dVdRa) - dZc * dVdRc;

    k.f1X = k.f1X - dXb * (-dVdRb) - dXc * (-dVdRc);
    k.f1Y = k.f1Y - dYb * (-dVdRb) - dYc * (-dVdRc);
    k.f1Z = k.f1Z - dZb * (-dVdRb) - dZc * (-dVdRc);
}