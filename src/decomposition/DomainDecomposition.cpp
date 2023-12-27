#include "DomainDecomposition.hpp"

DomainDecomposition::DomainDecomposition() {}
DomainDecomposition::~DomainDecomposition() {}

int DomainDecomposition::GetNumOfMyParticles() { return this->myParticles.size(); }
std::vector<Utility::Particle> DomainDecomposition::GetMyParticles() { return this->myParticles; }
void DomainDecomposition::SetMyParticles(std::vector<Utility::Particle>& particles)
{
    // TODO: avoid this copy
    this->myParticles = particles;
}
void DomainDecomposition::Init(std::shared_ptr<Simulation> simulation) { this->simulation = simulation; }

void DomainDecomposition::updateMyParticles(double dt, Eigen::Vector3d gForce)
{
    // update all my particles
    for (size_t i = 0; i < myParticles.size(); i++) {
        if (!myParticles[i].isDummy) {
            myParticles[i].Update(dt, gForce);
        }
    }
}

void DomainDecomposition::updateMyParticlesPredictorStage(double dt)
{
    for (size_t i = 0; i < myParticles.size(); i++) {
        if (!myParticles[i].isDummy) {
            myParticles[i].UpdatePredictorStage(dt);
        }
    }
}

void DomainDecomposition::ResetForces()
{
    for (size_t i = 0; i < myParticles.size(); i++) {
        myParticles[i].ResetForce();
    }
}