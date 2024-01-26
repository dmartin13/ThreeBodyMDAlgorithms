#include "Simulation.hpp"

Simulation::Simulation(Utility::cliArguments& args, int iterations, int respaStepSize,
                       std::shared_ptr<Algorithm> algorithm, std::shared_ptr<Topology> topology,
                       std::shared_ptr<PairwisePotential> pairwisepotential,
                       std::shared_ptr<TriwisePotential> triwisepotential,
                       std::shared_ptr<DomainDecomposition> decomposition, MPI_Datatype* mpiParticleType,
                       std::vector<Utility::Particle>& particles, double dt, Eigen::Vector3d gForce,
                       std::string csvOutput)
    : args(args), iterations(iterations), respaStepSize(respaStepSize), algorithm(algorithm), topology(topology),
      pairwisepotential(pairwisepotential), triwisepotential(triwisepotential), decomposition(decomposition),
      mpiParticleType(mpiParticleType), particles(particles), dt(dt), gForce(gForce), csvOutput(csvOutput) {}

void Simulation::Init() {
    // if add brownian motion is true we add it here to all particles. This is a really easy way so that all particles
    // always have the same random velocity regardless of num processors
    if (args.useThermostat and dt != 0) {
        if (args.addBrownianMotion) {
            Thermostat::addBrownianMotion(particles, args.targetTemperature);
        }
    }

    if (args.respaStepSize > 0) {
        if (args.iterations % args.respaStepSize != 0) {
            throw std::runtime_error("The number of iterations must be divisible by the respa step-size");
        }

        if (args.useThermostat) {
            if (args.thermostatInterval % args.respaStepSize != 0) {
                throw std::runtime_error("The thermostat interval must be divisible by the respa step-size");
            }
        }
    }

    std::shared_ptr<Simulation> simulationPtr = shared_from_this();
    this->topology->Init(simulationPtr);
    this->decomposition->Init(simulationPtr);
    this->algorithm->Init(simulationPtr);
    this->pairwisepotential->Init(simulationPtr);
    this->triwisepotential->Init(simulationPtr);
}

Simulation::~Simulation() {}

std::shared_ptr<Algorithm> Simulation::GetAlgorithm() { return this->algorithm; }
std::shared_ptr<Topology> Simulation::GetTopology() { return this->topology; }
std::shared_ptr<PairwisePotential> Simulation::GetPairwisePotential() { return this->pairwisepotential; }
std::shared_ptr<TriwisePotential> Simulation::GetTriwisePotential() { return this->triwisepotential; }
std::shared_ptr<DomainDecomposition> Simulation::GetDecomposition() { return this->decomposition; }

void Simulation::Start() {
    if (args.useThermostat and dt != 0) {
        // Set the simulation directly to the desired initial temperature.
        Thermostat::apply(decomposition->GetMyParticles(), args.initialTemperature, std::numeric_limits<double>::max());
    }

    const bool respaActive = respaStepSize > 0;

    for (int i = 0; i < iterations; ++i) {
        // write simulation step
#if !defined(BENCHMARK_3BMDA) && !defined(TESTS_3BMDA) && !defined(PROFILE_3BMDA)
        if (this->csvOutput.compare("") != 0 and (i % this->args.csvWriteInterval == 0)) {
            writeSimulationStepToCSV(this->csvOutput.substr(0, this->csvOutput.find_last_of('.')) + "_" +
                                     std::to_string(i) + ".csv");
        }
#endif

#ifdef MEASURESIMSTEP_3BMDA
        std::chrono::time_point<std::chrono::system_clock> start;
        std::chrono::time_point<std::chrono::system_clock> end;
        start = std::chrono::system_clock::now();
#endif
        // determine if this iteration is a respa step
        const bool isRespaIteration = respaActive and (i % respaStepSize == 0);
        const bool nextIsRespaIteration = respaActive and ((i + 1) % respaStepSize == 0);

        // if this is the first iteration we need a force calculation here so we can update the velocities with the
        // threebody force if respa is used. If respa is not used, this forces are used in the first position update of
        // the inner loop.
        if (i == 0) {
            this->algorithm->SimulationStep(ForceType::TwoAndThreeBody);
            MPI_Barrier(this->topology->GetComm());
            // we dont use the upot calculated in the initial force computation. We calculate the potential energy after
            // the first position update
            pairwisepotential->GetAndResetPotentialEnergy();
            triwisepotential->GetAndResetPotentialEnergy();
            pairwisepotential->GetAndResetVirial();
            triwisepotential->GetAndResetVirial();
        }

        if (respaActive and isRespaIteration) {
            // update velocities with three body force
            decomposition->UpdateVelocities(dt, ForceType::ThreeBody, respaStepSize);
        }

        // the following part is the inner loop of the respa algorithm. If respa is not used, this corresponds to the
        // störmer verlet integration
        {
            // update particle positions using the two-body force
            decomposition->UpdatePositions(dt, gForce, respaActive ? ForceType::TwoBody : ForceType::TwoAndThreeBody);
            MPI_Barrier(this->topology->GetComm());

            // execute algorithm... force calculation
            // determine the forceType to calculate in this iteration
            auto forceTypeToCalculate = ForceType::TwoAndThreeBody;
            if (respaActive and (not nextIsRespaIteration)) {
                // if respa should be used but the next iteration is not a respa iteration then only calculate two body
                // interactions
                forceTypeToCalculate = ForceType::TwoBody;
            }

            this->algorithm->SimulationStep(forceTypeToCalculate);
            MPI_Barrier(this->topology->GetComm());

            // update the velocities
            decomposition->UpdateVelocities(dt, respaActive ? ForceType::TwoBody : ForceType::TwoAndThreeBody,
                                            respaStepSize);

            if (not respaActive) {
                // apply thermostat
                if (args.useThermostat) {
                    if ((i + 1) % (args.thermostatInterval) == 0) {
                        Thermostat::apply(decomposition->GetMyParticles(), args.targetTemperature,
                                          args.deltaTemperature);
                        // std::cout << "apply thermostat in iteration " << i << std::endl;
                    }
                }

                // calculate kinetic energy
                const auto eKin = calculateKineticEnergy();
                kineticEnergy.push_back(eKin);

                // get the potential energy
                const auto uPotTwoBody = pairwisepotential->GetAndResetPotentialEnergy();
                const auto uPotThreeBody = triwisepotential->GetAndResetPotentialEnergy();
                potentialEnergyTwoBody.push_back(uPotTwoBody);
                potentialEnergyThreeBody.push_back(uPotThreeBody);

                // total energy
                totalEnergy.push_back(uPotTwoBody + uPotThreeBody + eKin);

            } else if (not nextIsRespaIteration) {
                pairwisepotential->GetAndResetPotentialEnergy();
                triwisepotential->GetAndResetPotentialEnergy();
            }

            MPI_Barrier(this->topology->GetComm());
        }

        // check if the next iteration is a respa iteration
        if (respaActive and nextIsRespaIteration) {
            // 3-body force has already bee calculated in inner loop
            // update velocities using three body force
            decomposition->UpdateVelocities(dt, ForceType::ThreeBody, respaStepSize);

            // apply thermostat
            if (args.useThermostat) {
                if ((i + 1) % (args.thermostatInterval) == 0) {
                    Thermostat::apply(decomposition->GetMyParticles(), args.targetTemperature, args.deltaTemperature);
                    // std::cout << "apply thermostat in iteration " << i << std::endl;
                }
            }

            // calculate kinetic energy
            const auto eKin = calculateKineticEnergy();
            kineticEnergy.push_back(eKin);

            // get the potential energy
            const auto uPotTwoBody = pairwisepotential->GetAndResetPotentialEnergy();
            const auto uPotThreeBody = triwisepotential->GetAndResetPotentialEnergy();
            potentialEnergyTwoBody.push_back(uPotTwoBody);
            potentialEnergyThreeBody.push_back(uPotThreeBody);

            // total energy
            totalEnergy.push_back(uPotTwoBody + uPotThreeBody + eKin);
        }

        // reset virial
        pairwisepotential->GetAndResetVirial();
        triwisepotential->GetAndResetVirial();

#ifdef MEASURESIMSTEP_3BMDA
        end = std::chrono::system_clock::now();
        auto elapsed_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
        std::vector<uint64_t> allTimes;
        if (topology->GetWorldRank() == 0) {
            allTimes.resize(topology->GetWorldSize());
        }
        MPI_Gather(&elapsed_time, 1, MPI_INT64_T, allTimes.data(), 1, MPI_INT64_T, 0, topology->GetComm());
        if (topology->GetWorldRank() == 0) {
            uint64_t max = *max_element(allTimes.begin(), allTimes.end());
            std::string allTimesStr = "[";
            for (size_t i = 0; i < allTimes.size(); i++) {
                allTimesStr.append(std::to_string(allTimes[i]));
                if (i < allTimes.size() - 1) {
                    allTimesStr.append(", ");
                }
            }
            allTimesStr.append("]");
            std::string json = "{\"allTimes\": " + allTimesStr + ", \"max\": " + std::to_string(max) + "}";
            std::cout << json << std::endl;
        }

#endif

        // print out some statistics of the simulation
        if ((not respaActive) or (respaActive and nextIsRespaIteration)) {
            double potentialEnergyAfterIterationTwoBody = potentialEnergyTwoBody.back();
            double potentialEnergyAfterIterationThreeBody = potentialEnergyThreeBody.back();
            double kineticEnergyAfterIteration = kineticEnergy.back();
            if (this->GetTopology()->GetWorldRank() == 0) {
                MPI_Reduce(MPI_IN_PLACE, &potentialEnergyAfterIterationTwoBody, 1, MPI_DOUBLE, MPI_SUM, 0,
                           topology->GetComm());
                MPI_Reduce(MPI_IN_PLACE, &potentialEnergyAfterIterationThreeBody, 1, MPI_DOUBLE, MPI_SUM, 0,
                           topology->GetComm());
                MPI_Reduce(MPI_IN_PLACE, &kineticEnergyAfterIteration, 1, MPI_DOUBLE, MPI_SUM, 0, topology->GetComm());
                std::cout << "Finalized iteration         : " << i << "\n"
                          << "Potential Energy Two-Body   : " << std::format("{}", potentialEnergyAfterIterationTwoBody)
                          << "\n"
                          << "Potential Energy Three-Body : "
                          << std::format("{}", potentialEnergyAfterIterationThreeBody) << "\n"
                          << "Kinetic Energy              : " << std::format("{}", kineticEnergyAfterIteration)
                          << std::endl;
            } else {
                MPI_Reduce(&potentialEnergyAfterIterationTwoBody, nullptr, 1, MPI_DOUBLE, MPI_SUM, 0,
                           topology->GetComm());
                MPI_Reduce(&potentialEnergyAfterIterationThreeBody, nullptr, 1, MPI_DOUBLE, MPI_SUM, 0,
                           topology->GetComm());
                MPI_Reduce(&kineticEnergyAfterIteration, nullptr, 1, MPI_DOUBLE, MPI_SUM, 0, topology->GetComm());
            }
        }
    }

    // Record last state of simulation.
#if !defined(BENCHMARK_3BMDA) && !defined(TESTS_3BMDA) && !defined(PROFILE_3BMDA)
    if (this->csvOutput.compare("") != 0) {
        writeSimulationStepToCSV(this->csvOutput.substr(0, this->csvOutput.find_last_of('.')) + "_" +
                                 std::to_string(iterations) + ".csv");
    }
#endif
}

double Simulation::calculateKineticEnergy() {
    double kinEAcc = 0;
    for (const auto& p : decomposition->GetMyParticles()) {
        if (p.isDummy) {
            continue;
        }
        kinEAcc += p.mass * (p.vX * p.vX + p.vY * p.vY + p.vZ * p.vZ);
    }
    kinEAcc *= 0.5;
    return kinEAcc;
}

MPI_Datatype* Simulation::GetMPIParticleType() { return this->mpiParticleType; }

std::vector<Utility::Particle>& Simulation::GetAllParticles() { return this->particles; }

double Simulation::GetDeltaT() { return this->dt; }

int Simulation::GetNumIterations() { return this->iterations; }

uint64_t Simulation::GetNumBufferInteractions(int step) {
    return (size_t)step < this->numInteractions.size() ? std::get<0>(this->numInteractions[step]) : 0;
}

uint64_t Simulation::GetNumParticleInteractions(int step) {
    return (size_t)step < this->numInteractions.size() ? std::get<1>(this->numInteractions[step]) : 0;
}

Eigen::Vector3d Simulation::GetGForce() { return this->gForce; }

std::vector<double>& Simulation::GetKineticEnergy() { return kineticEnergy; }
std::vector<double>& Simulation::GetPotentialEnergyTwoBody() { return potentialEnergyTwoBody; }
std::vector<double>& Simulation::GetPotentialEnergyThreeBody() { return potentialEnergyThreeBody; }
std::vector<double>& Simulation::GetTotalEnergy() { return totalEnergy; }

void Simulation::writeSimulationStepToCSV(std::string file) {
    std::vector<Utility::Particle> receivedParticlesForCSVOutput;

    int numProcessors = topology->GetWorldSize();
    std::vector<int> numParticlesPerProcessor;
    numParticlesPerProcessor.resize(numProcessors);
    // elements from each process are gathered in order of their rank
    int numOfMyParticles = decomposition->GetNumOfMyParticles();
    MPI_Gather(&numOfMyParticles, 1, MPI_INT, numParticlesPerProcessor.data(), 1, MPI_INT, 0, topology->GetComm());

    std::vector<int> displacements;
    int sumDispl = 0;
    for (int i = 0; i < numProcessors; i++) {
        displacements.push_back(sumDispl);
        sumDispl += numParticlesPerProcessor[i];
    }

    if (topology->GetWorldRank() == 0) {
        receivedParticlesForCSVOutput.resize(sumDispl);
    }

    std::vector<Utility::Particle> particlesToSend = this->decomposition->GetMyParticles();

    MPI_Gatherv(particlesToSend.data(), particlesToSend.size(), *mpiParticleType, receivedParticlesForCSVOutput.data(),
                numParticlesPerProcessor.data(), displacements.data(), *mpiParticleType, 0, topology->GetComm());

    std::sort(receivedParticlesForCSVOutput.begin(), receivedParticlesForCSVOutput.end(),
              [](Utility::Particle a, Utility::Particle b) { return a.ID < b.ID; });

    if (topology->GetWorldRank() == 0) {
        Utility::writeStepToCSVWithForces(file, receivedParticlesForCSVOutput);
    }
}