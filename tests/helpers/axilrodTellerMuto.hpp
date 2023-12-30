#pragma once
#include <Eigen/Core>
#include <tuple>

std::tuple<Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d> calculateATMForce(const Eigen::Vector3d& posI,
                                                                                const Eigen::Vector3d& posJ,
                                                                                const Eigen::Vector3d& posK,
                                                                                const double nu) {
    const auto displacementIJ = posJ - posI;
    const auto displacementJK = posK - posJ;
    const auto displacementKI = posI - posK;

    const double distSquaredIJ = displacementIJ.dot(displacementIJ);
    const double distSquaredJK = displacementJK.dot(displacementJK);
    const double distSquaredKI = displacementKI.dot(displacementKI);

    // Dot products of both distance vectors going from one particle
    const double IJDotKI = displacementIJ.dot(displacementKI);
    const double IJDotJK = displacementIJ.dot(displacementJK);
    const double JKDotKI = displacementJK.dot(displacementKI);

    const double allDotProducts = IJDotKI * IJDotJK * JKDotKI;

    const double allDistsSquared = distSquaredIJ * distSquaredJK * distSquaredKI;
    const double allDistsTo5 = allDistsSquared * allDistsSquared * std::sqrt(allDistsSquared);
    const double factor = 3.0 * nu / allDistsTo5;

    const auto forceIDirectionJK = displacementJK * IJDotKI * (IJDotJK - JKDotKI);
    const auto forceIDirectionIJ =
        displacementIJ * (IJDotJK * JKDotKI - distSquaredJK * distSquaredKI + 5.0 * allDotProducts / distSquaredIJ);
    const auto forceIDirectionKI =
        displacementKI * (-IJDotJK * JKDotKI + distSquaredIJ * distSquaredJK - 5.0 * allDotProducts / distSquaredKI);

    const auto forceI = (forceIDirectionJK + forceIDirectionIJ + forceIDirectionKI) * factor;

    const auto forceJDirectionKI = displacementKI * IJDotJK * (JKDotKI - IJDotKI);
    const auto forceJDirectionIJ =
        displacementIJ * (-IJDotKI * JKDotKI + distSquaredJK * distSquaredKI - 5.0 * allDotProducts / distSquaredIJ);
    const auto forceJDirectionJK =
        displacementJK * (IJDotKI * JKDotKI - distSquaredIJ * distSquaredKI + 5.0 * allDotProducts / distSquaredJK);

    const auto forceJ = (forceJDirectionKI + forceJDirectionIJ + forceJDirectionJK) * factor;

    const auto forceK = (forceI + forceJ) * (-1.0);

    return {forceI, forceJ, forceK};
}