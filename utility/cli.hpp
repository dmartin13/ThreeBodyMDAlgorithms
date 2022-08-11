#include <iostream>
#include <string>
#include <vector>

#include "enums.hpp"

#pragma once

namespace Utility
{
    struct cliArguments {
        int iterations;
        double cutoff;
        double deltaT;
        std::string inputCSV;
        AlgorithmType algorithm;
        Eigen::Vector3d gForce;

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
                        a.printHelp();
                        exit(1);
                    }
                } else {
                    a.printHelp();
                    exit(1);
                }
            } catch (std::exception& e) {
                a.printHelp();
                exit(1);
            }
        }
        return a;
    }
}  // namespace Utility