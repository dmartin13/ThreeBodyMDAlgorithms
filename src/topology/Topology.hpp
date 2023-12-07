#pragma once

#include <mpi.h>

#include <memory>

#include "fwd.hpp"
#include "simulation/Simulation.hpp"
#include "utility/utility.hpp"

class Topology {
protected:
    int worldRank;
    int worldSize;
    std::shared_ptr<Simulation> simulation;
    MPI_Comm comm;

public:
    Topology();
    virtual ~Topology();
    int GetWorldRank();
    int GetWorldSize();
    virtual void Init(std::shared_ptr<Simulation> simulation);
    MPI_Comm GetComm();
};