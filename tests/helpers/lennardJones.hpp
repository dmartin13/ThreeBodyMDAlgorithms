#pragma once
#include <Eigen/Core>
#include <tuple>

std::tuple<Eigen::Vector3d, Eigen::Vector3d> calculateLJForce(const Eigen::Vector3d &posI, const Eigen::Vector3d &posJ,
                                                              const double sigma, const double epsilon) {
    // the vector from particle j to i
    const auto posIMinusPosJ = posI - posJ;

    // the distance between both particles
    const auto dist = std::sqrt(posIMinusPosJ.dot(posIMinusPosJ));

    // r^6
    const auto r6 = dist * dist * dist * dist * dist * dist;
    // sigma^6
    const auto sigma6 = sigma * sigma * sigma * sigma * sigma * sigma;
    // sigma^6 / r^7
    const auto dlj6 = sigma6 / (r6 * dist);
    // sigma^12 / r^13
    const auto dlj12 = (sigma6 * sigma6) / (r6 * r6 * dist);
    // the derivative with respect to r of the lennard jones potential
    const auto dUr = 48. * epsilon * (dlj12 - 0.5 * dlj6);

    // the forces in x, y and z direction
    const auto fx = (posIMinusPosJ[0] / dist) * dUr;
    const auto fy = (posIMinusPosJ[1] / dist) * dUr;
    const auto fz = (posIMinusPosJ[2] / dist) * dUr;

    const auto forceI = Eigen::Vector3d{fx, fy, fz};
    const auto forceJ = Eigen::Vector3d{-fx, -fy, -fz};

    return {forceI, forceJ};
}