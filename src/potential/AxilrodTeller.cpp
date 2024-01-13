#include "AxilrodTeller.hpp"

AxilrodTeller::AxilrodTeller(double nu) : nu(nu) {}

void AxilrodTeller::CalculateForces(Utility::Particle &i, Utility::Particle &j, Utility::Particle &k) {
    // we assume forces are set to 0 for each particle

    const Eigen::Vector3d displacementIJ = j.GetR() - i.GetR();
    const Eigen::Vector3d displacementJK = k.GetR() - j.GetR();
    const Eigen::Vector3d displacementKI = i.GetR() - k.GetR();

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

    i.f1X += forceI[0];
    i.f1Y += forceI[1];
    i.f1Z += forceI[2];

    j.f1X += forceJ[0];
    j.f1Y += forceJ[1];
    j.f1Z += forceJ[2];

    k.f1X += forceK[0];
    k.f1Y += forceK[1];
    k.f1Z += forceK[2];

    // potential energy for this triple interaction
    const double upot = factor * (allDistsSquared - 3.0 * allDotProducts) / 3.0;
    potentialEnergy += upot;

    // Virial is calculated as f_i * r_i
    // see Thompson et al.: https://doi.org/10.1063/1.3245303
    virial += (Eigen::Array3d{forceI[0], forceI[1], forceI[2]} * i.GetR() +
               Eigen::Array3d{forceJ[0], forceJ[1], forceJ[2]} * j.GetR() +
               Eigen::Array3d{forceK[0], forceK[1], forceK[2]} * k.GetR());

    perfomedInteractions++;
}