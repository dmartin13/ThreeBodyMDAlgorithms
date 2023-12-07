#pragma once

#include <memory>

#include "../fwd.hpp"
#include "../utility/utility.hpp"
#include "Potential.hpp"

class PairwisePotential : public Potential {
protected:
public:
    virtual void CalculateForces(Utility::Particle &i, Utility::Particle &j) = 0;
};