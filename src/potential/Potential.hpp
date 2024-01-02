#pragma once

#include <memory>

#include "fwd.hpp"
#include "utility/utility.hpp"

#ifdef PROFILE_3BMDA
#include <chrono>
#endif

class Potential {
protected:
    std::shared_ptr<Simulation> simulation;
    double potentialEnergy;
    Eigen::Array3d virial;

public:
    Potential();
    virtual void Init(std::shared_ptr<Simulation> simulation);
    double GetAndResetPotentialEnergy();
    Eigen::Array3d GetAndResetVirial();
};