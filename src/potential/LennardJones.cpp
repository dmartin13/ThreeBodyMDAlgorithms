#include "LennardJones.hpp"

LennardJones::LennardJones(double epsilon, double sigma)
    : _epsilon(epsilon), _sigma(sigma), _sigmaSquared(sigma * sigma), _epsilon24(epsilon * 24) {}

void LennardJones::CalculateForces(Utility::Particle &i, Utility::Particle &j) {
    auto sigmaSquared = _sigmaSquared;
    auto epsilon24 = _epsilon24;
    Eigen::Array3d dr = i.GetR() - j.GetR();
    double dr2 = dr.cwiseProduct(dr).sum();

    double invdr2 = 1. / dr2;
    double lj6 = sigmaSquared * invdr2;
    lj6 = lj6 * lj6 * lj6;
    double lj12 = lj6 * lj6;
    double lj12m6 = lj12 - lj6;
    double fac = epsilon24 * (lj12 + lj12m6) * invdr2;
    auto f = dr * fac;

    i.f0X += f.x();
    i.f0Y += f.y();
    i.f0Z += f.z();

    // we use newton 3
    j.f0X -= f.x();
    j.f0Y -= f.y();
    j.f0Z -= f.z();

    // calculate potential energy and virial
    potentialEnergy += 4.0 * _epsilon * lj12m6;
    virial += dr * f;

    perfomedInteractions++;
}