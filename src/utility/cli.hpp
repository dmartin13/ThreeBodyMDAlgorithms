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
        double cutoff;
        double deltaT;
        std::string inputCSV;
        std::string outputCSV;
        AlgorithmType algorithm;
        Eigen::Vector3d gForce;
        bool optimalDecomposition;
        std::string outputProfile;
        std::string benchYaml;

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