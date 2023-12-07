#pragma once

#include <memory>

#include "../fwd.hpp"
#include "../utility/utility.hpp"
#include "Potential.hpp"

class TriwisePotential : public Potential {
protected:
public:
    virtual void CalculateForces(Utility::Particle &i, Utility::Particle &j, Utility::Particle &k) = 0;
};