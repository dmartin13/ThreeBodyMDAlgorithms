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

#ifdef PROFILE_3BMDA
    uint64_t timeAcc;
    size_t counter;
#endif

public:
    virtual void Init(std::shared_ptr<Simulation> simulation);

#ifdef PROFILE_3BMDA
    std::map<std::string, std::pair<char, double>> GetAvgCalcTime();
    void ResetTime();
#endif
};