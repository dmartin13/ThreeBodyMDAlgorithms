#include "DomainDecomposition.hpp"

DomainDecomposition::DomainDecomposition() {}
DomainDecomposition::~DomainDecomposition() {}

int DomainDecomposition::GetNumOfMyParticles() { return this->myParticles.size(); }
std::vector<Utility::Particle>& DomainDecomposition::GetMyParticles() { return this->myParticles; }
void DomainDecomposition::SetMyParticles(std::vector<Utility::Particle>& particles) {
    // TODO: avoid this copy
    this->myParticles = particles;
}
void DomainDecomposition::Init(std::shared_ptr<Simulation> simulation) { this->simulation = simulation; }

void DomainDecomposition::UpdatePositions(double dt, Eigen::Vector3d gForce, ForceType forceTypeToUse) {
    TimeIntegration::calculatePositionsAndResetForces(myParticles, dt, {gForce[0], gForce[1], gForce[2]},
                                                      forceTypeToUse);
}

void DomainDecomposition::UpdateVelocities(double dt, ForceType forceType, size_t respaStepSize) {
    if (dt == 0) {
        return;
    }

    if (forceType == ForceType::ThreeBody) {
        TimeIntegration::calculateVelocities(myParticles, dt, /*outerRespaStep*/ true, respaStepSize);
    } else if (forceType == ForceType::TwoBody) {
        TimeIntegration::calculateVelocities(myParticles, dt, /*outerRespaStep*/ false, respaStepSize);
    } else {
        TimeIntegration::calculateVelocities(myParticles, dt, /*outerRespaStep*/ false, -1);
    }
}