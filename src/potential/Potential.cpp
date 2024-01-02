#include "Potential.hpp"

Potential::Potential() : potentialEnergy(0), virial(0) {}

double Potential::GetAndResetPotentialEnergy() {
    const auto retVal = potentialEnergy;
    potentialEnergy = 0;
    return retVal;
}

Eigen::Array3d Potential::GetAndResetVirial() {
    const auto retVal = virial;
    virial = {0, 0, 0};
    return retVal;
}

void Potential::Init(std::shared_ptr<Simulation> simulation) {
    this->simulation = simulation;
}