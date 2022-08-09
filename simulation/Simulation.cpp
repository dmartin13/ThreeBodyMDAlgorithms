#include "Simulation.hpp"

Simulation::Simulation(int iterations, std::shared_ptr<Algorithm> algorithm, std::shared_ptr<Topology> topology,
                       std::shared_ptr<Potential> potential, std::shared_ptr<DomainDecomposition> decomposition,
                       MPI_Datatype* mpiParticleType, std::vector<Utility::Particle>& particles)
    : iterations(iterations), algorithm(algorithm), topology(topology), potential(potential),
      decomposition(decomposition), mpiParticleType(mpiParticleType), particles(particles)
{}

void Simulation::Init()
{
    std::shared_ptr<Simulation> simulationPtr = shared_from_this();
    this->topology->Init(simulationPtr);
    this->decomposition->Init(simulationPtr);
    this->algorithm->Init(simulationPtr);
    this->potential->Init(simulationPtr);
}

Simulation::~Simulation() {}

std::shared_ptr<Algorithm> Simulation::GetAlgorithm() { return this->algorithm; }
std::shared_ptr<Topology> Simulation::GetTopology() { return this->topology; }
std::shared_ptr<Potential> Simulation::GetPotential() { return this->potential; }
std::shared_ptr<DomainDecomposition> Simulation::GetDecomposition() { return this->decomposition; }

void Simulation::Start()
{
    for (int i = 0; i < iterations; ++i) {
        // do step
        this->algorithm->SimulationStep();
        //this->decomposition->Update();
    }
}

MPI_Datatype* Simulation::GetMPIParticleType() { return this->mpiParticleType; }

std::vector<Utility::Particle>& Simulation::GetAllParticles() { return this->particles; }