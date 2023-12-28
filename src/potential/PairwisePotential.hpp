#pragma once

#include <memory>

#include "Potential.hpp"
#include "fwd.hpp"
#include "utility/utility.hpp"

class PairwisePotential : public Potential {
protected:
public:
    virtual void CalculateForces(Utility::Particle &i, Utility::Particle &j) = 0;
};