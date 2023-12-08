#pragma once

#include <Eigen/Dense>
#include <iostream>
#include <string>
#include <vector>

#include "enums.hpp"

namespace Utility
{
    struct cliArguments {
        int iterations;
        double cutoff;
        double deltaT;
        std::string inputCSV;
        std::string outputCSV;
        AlgorithmType algorithm;
        Eigen::Vector3d gForce;
        bool optimalDecomposition;
        std::string outputProfile;
        std::string benchYaml;

        void printHelp();
    };

    void cliArguments::printHelp()
    {
        std::cout << "Help: \n"
                  << "Options:\n"
                  << "\t-h,--help\t\tShow this help message\n"
                  << "\t-a,--algorithm\t\talgorithm to use (\"nata\", \"p3bca\", \"auta\")\n"
                  << "\t-i,--iterations\t\tnum of iterations to simulate\n"
                  << "\t-d,--delta\t\tduration of one simulation step\n"
                  << "\t-gx,--gravityZ\t\tgravitational force in x-direction\n"
                  << "\t-gy,--gravityY\t\tgravitational force in y-direction\n"
                  << "\t-gz,--gravityZ\t\tgravitational force in z-direction\n"
                  << "\t-csv,--csv\t\tcsv file with particles\n"
                  << "\t-o,--out\t\tcsv base file name that will be used to create all csv outputs from each "
                     "simulation step\n"
                  << "-op,--outProfile\t\tjson file to save profile output\n"
                  << "-decomp,--decomp\t\t(optimal | naive) decomposition strategy for cart topology\n"
                  << "-benchyaml,--benchyaml\t\tYAML file that contains config for benchmarks\n"
                  << "\t-c,--cutoff\t\tcutoff distance" << std::endl;
    }

    cliArguments cliParse(std::vector<std::string> args)
    {
        cliArguments a;
        std::string flag, value;

        for (int i = 0; (size_t)i < args.size(); i++) {
            try {
                if (args[i].rfind("-", 0) == 0) {
                    flag = args[i].substr(1);
                    if (flag.compare("c") == 0 || flag.compare("-cutoff") == 0) {
                        if (args.size() <= (size_t)(i + 1)) {
                            a.printHelp();
                            exit(1);
                        }
                        value = args[i + 1];
                        a.cutoff = std::stod(value);
                        i++;
                    } else if (flag.compare("i") == 0 || flag.compare("-iterations") == 0) {
                        if (args.size() <= (size_t)(i + 1)) {
                            a.printHelp();
                            exit(1);
                        }
                        value = args[i + 1];
                        a.iterations = std::stoi(value);
                        i++;
                    } else if (flag.compare("a") == 0 || flag.compare("-algorithm") == 0) {
                        if (args.size() <= (size_t)(i + 1)) {
                            a.printHelp();
                            exit(1);
                        }
                        value = args[i + 1];
                        if (value.compare("nata") == 0) {
                            a.algorithm = AlgorithmType::NATAType;
                        } else if (value.compare("p3bca") == 0) {
                            a.algorithm = AlgorithmType::P3BCAType;
                        } else if (value.compare("auta") == 0) {
                            a.algorithm = AlgorithmType::AUTAType;
                        } else {
                            a.printHelp();
                            exit(1);
                        }
                        i++;
                    } else if (flag.compare("decomp") == 0 || flag.compare("-decomp") == 0) {
                        if (args.size() <= (size_t)(i + 1)) {
                            a.printHelp();
                            exit(1);
                        }
                        value = args[i + 1];
                        if (value.compare("optimal") == 0) {
                            a.optimalDecomposition = true;
                        } else if (value.compare("naive") == 0) {
                            a.optimalDecomposition = false;
                        } else {
                            a.printHelp();
                            exit(1);
                        }
                        i++;
                    } else if (flag.compare("d") == 0 || flag.compare("-delta") == 0) {
                        if (args.size() <= (size_t)(i + 1)) {
                            a.printHelp();
                            exit(1);
                        }
                        value = args[i + 1];
                        a.deltaT = std::stod(value);
                        i++;
                    } else if (flag.compare("csv") == 0 || flag.compare("-csv") == 0) {
                        if (args.size() <= (size_t)(i + 1)) {
                            a.printHelp();
                            exit(1);
                        }
                        value = args[i + 1];
                        a.inputCSV = value;
                        i++;
                    } else if (flag.compare("o") == 0 || flag.compare("-out") == 0) {
                        if (args.size() <= (size_t)(i + 1)) {
                            a.printHelp();
                            exit(1);
                        }
                        value = args[i + 1];
                        a.outputCSV = value;
                        i++;
                    } else if (flag.compare("op") == 0 || flag.compare("-outProfile") == 0) {
                        if (args.size() <= (size_t)(i + 1)) {
                            a.printHelp();
                            exit(1);
                        }
                        value = args[i + 1];
                        a.outputProfile = value;
                        i++;
                    } else if (flag.compare("benchyaml") == 0 || flag.compare("-benchyaml") == 0) {
                        if (args.size() <= (size_t)(i + 1)) {
                            a.printHelp();
                            exit(1);
                        }
                        value = args[i + 1];
                        a.benchYaml = value;
                        i++;
                    } else if (flag.compare("gx") == 0 || flag.compare("-gravityX") == 0) {
                        if (args.size() <= (size_t)(i + 1)) {
                            a.printHelp();
                            exit(1);
                        }
                        value = args[i + 1];
                        a.gForce[0] = std::stod(value);

                        i += 1;
                    } else if (flag.compare("gy") == 0 || flag.compare("-gravityY") == 0) {
                        if (args.size() <= (size_t)(i + 1)) {
                            a.printHelp();
                            exit(1);
                        }
                        value = args[i + 1];
                        a.gForce[1] = std::stod(value);

                        i += 1;
                    } else if (flag.compare("gz") == 0 || flag.compare("-gravityZ") == 0) {
                        if (args.size() <= (size_t)(i + 1)) {
                            a.printHelp();
                            exit(1);
                        }
                        value = args[i + 1];
                        a.gForce[2] = std::stod(value);

                        i += 1;
                    } else {
                        std::cout << "CLI parser: unrecognized option: " << flag << std::endl;
                        // a.printHelp();
                        // exit(1);
                    }
                } else {
                    std::cout << "CLI parser: unrecognized option: " << args[i] << std::endl;
                    // a.printHelp();
                    // exit(1);
                }
            } catch (std::exception& e) {
                a.printHelp();
                exit(1);
            }
        }
        return a;
    }
}  // namespace Utility