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
        bool useThermostat;
        bool addBrownianMotion = false;
        size_t thermostatInterval;
        double initialTemperature;
        double targetTemperature;
        double deltaTemperature;

        void printHelp();
    };

    cliArguments cliParse(std::vector<std::string> args);
}  // namespace Utility