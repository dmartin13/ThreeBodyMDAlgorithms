#pragma once

#include <Eigen/Dense>
#include <iostream>
#include <string>
#include <vector>

#include "enums.hpp"

namespace Utility {
    struct cliArguments {
        int iterations;
        int respaStepSize;
        size_t csvWriteInterval{1};
        double cutoff;
        double deltaT;
        std::string inputCSV;
        std::string outputCSV;
        AlgorithmType algorithm;
        Eigen::Vector3d gForce;
        bool optimalDecomposition;
        std::string outputProfile;
        std::string benchYaml;
        unsigned reflect;
        bool disableThreebodyInteractions{false};
        std::array<double, 3> boxSize;
        std::array<double, 3> bottomLeft{std::numeric_limits<double>::infinity(),
                                         std::numeric_limits<double>::infinity(),
                                         std::numeric_limits<double>::infinity()};

        // thermostat parameters
        bool useThermostat{false};
        bool addBrownianMotion{false};
        size_t thermostatInterval{0};
        double initialTemperature{0.0};
        double targetTemperature{0.0};
        double deltaTemperature{0.0};

        // parameters for potentials
        double epsilon{1.0};
        double sigma{1.0};
        double nu{1.0};

        void printHelp();
    };

    cliArguments cliParse(std::vector<std::string> args);
}  // namespace Utility