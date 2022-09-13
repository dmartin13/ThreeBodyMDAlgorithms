#include "Simulation.hpp"

Simulation::Simulation(int iterations, std::shared_ptr<Algorithm> algorithm, std::shared_ptr<Topology> topology,
                       std::shared_ptr<Potential> potential, std::shared_ptr<DomainDecomposition> decomposition,
                       MPI_Datatype* mpiParticleType, std::vector<Utility::Particle>& particles, double dt,
                       Eigen::Vector3d gForce)
    : iterations(iterations), algorithm(algorithm), topology(topology), potential(potential),
      decomposition(decomposition), mpiParticleType(mpiParticleType), particles(particles), dt(dt), gForce(gForce)
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
        // do step and record the number of interactions
        numInteractions.push_back(this->algorithm->SimulationStep());
        MPI_Barrier(this->topology->GetComm());
        this->decomposition->Update(this->dt, this->gForce);
        MPI_Barrier(this->topology->GetComm());
    }
}

MPI_Datatype* Simulation::GetMPIParticleType() { return this->mpiParticleType; }

std::vector<Utility::Particle>& Simulation::GetAllParticles() { return this->particles; }

double Simulation::GetDeltaT() { return this->dt; }
int Simulation::GetNumIterations() { return this->iterations; }
int Simulation::GetNumBufferInteractions(int step)
{
    return (size_t)step < this->numInteractions.size() ? std::get<0>(this->numInteractions[step]) : 0;
}
int Simulation::GetNumParticleInteractions(int step)
{
    return (size_t)step < this->numInteractions.size() ? std::get<1>(this->numInteractions[step]) : 0;
}
Eigen::Vector3d Simulation::GetGForce() { return this->gForce; }