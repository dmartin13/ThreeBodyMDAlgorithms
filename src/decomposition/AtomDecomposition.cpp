#include "AtomDecomposition.hpp"

AtomDecomposition::AtomDecomposition() {}
AtomDecomposition::~AtomDecomposition() {}

void AtomDecomposition::Init(std::shared_ptr<Simulation> simulation) {
    DomainDecomposition::Init(simulation);

    this->worldRank = this->simulation->GetTopology()->GetWorldRank();
    this->worldSize = this->simulation->GetTopology()->GetWorldSize();

    std::vector<Utility::Particle>& particles = this->simulation->GetAllParticles();

    int numParticles = (int)particles.size();
    int offset;

    int numProcNormal = worldSize - (numParticles % worldSize);
    int numParticlesProcNormal = numParticles / worldSize;
    int numParticlesLastProcs = numParticlesProcNormal + 1;

    int numOfMyParticles;

    if (worldRank < numProcNormal) {
        numOfMyParticles = numParticlesProcNormal;
        offset = worldRank * numParticlesProcNormal;
    } else {
        numOfMyParticles = numParticlesLastProcs;
        offset = numProcNormal * numParticlesProcNormal + (worldRank - numProcNormal) * numOfMyParticles;
    }

    myParticles.clear();

    for (int i = offset; i < offset + numOfMyParticles; ++i) {
        myParticles.push_back(particles[i]);
    }

    // add dummy particles for padding
    if (worldRank < numProcNormal && (numParticles % worldSize) != 0) {
        myParticles.push_back(Utility::Particle(true));
    }
}