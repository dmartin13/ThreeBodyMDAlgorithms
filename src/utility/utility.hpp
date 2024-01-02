#pragma once

#include <memory>

#include "../external/rapidcsv/src/rapidcsv.h"
#include "structs.hpp"

namespace Utility {
    // https://stackoverflow.com/a/4609795
    template <typename T>
    static int sgn(T val);

    int mod(int a, int b);

    void getParticlesFromCSV(std::string file, std::vector<Particle> &particles);
    void getParticlesFromTuple(
        std::vector<std::tuple<int, double, double, double, double, double, double, double>> &tuples,
        std::vector<Particle> &particles);

    void writeStepToCSV(std::string file, std::vector<Particle> &particles);
    void writeStepToCSVWithForces(std::string file, std::vector<Particle> &particles);

    int BinomialCoefficient(const int n, const int k);

    std::string get_file_contents(const char *filename);

}  // namespace Utility