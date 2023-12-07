#include "Potential.hpp"

void Potential::Init(std::shared_ptr<Simulation> simulation)
{
    this->simulation = simulation;
#ifdef PROFILE_3BMDA
    this->counter = 0;
    this->timeAcc = 0;
#endif
}

#ifdef PROFILE_3BMDA
std::map<std::string, std::pair<char, double>> Potential::GetAvgCalcTime()
{
    std::map<std::string, std::pair<char, double>> time;

    bool hasKey = time.count("CalculateForces");
    if (!hasKey) {
        time["CalculateForces"] = std::make_pair(0, 0);
    }
    time["CalculateForces"].second = (double)this->timeAcc / (double)this->counter;

    return time;
}
void Potential::ResetTime()
{
    this->counter = 0;
    this->timeAcc = 0;
}
#endif