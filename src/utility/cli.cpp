#include "cli.hpp"

namespace Utility {

    void cliArguments::printHelp() {
        std::cout << "Help: \n"
                  << "Options:\n"
                  << "\t-h,--help\t\tShow this help message\n"
                  << "\t-a,--algorithm\t\talgorithm to use (\"eauta\")\n"
                  << "\t-i,--iterations\t\tnum of iterations to simulate\n"
                  << "\t-r,--respaStepSize\tr-RESPA Step Size\n"
                  << "\t-tI,--thermostatInterval\tThermostat Interval\n"
                  << "\t-tInitT,--thermostatInitialTemperature\tThermostat Initial Temperature\n"
                  << "\t-tDeltaT,--thermostatDeltaTemperature\tThermostat Delta Temperature\n"
                  << "\t-tTargetT,--thermostatTargetTemperature\tThermostat Target Temperature\n"
                  << "\t-tAddB,--thermostatAddBrownianMotion\tThermostat Add Brownian Motion\n"
                  << "\t-d,--delta\t\tduration of one simulation step\n"
                  << "\t-gx,--gravityZ\t\tgravitational force in x-direction\n"
                  << "\t-gy,--gravityY\t\tgravitational force in y-direction\n"
                  << "\t-gz,--gravityZ\t\tgravitational force in z-direction\n"
                  << "\t-s,--sigma\t\tsigma for LJ potential\n"
                  << "\t-e,--epsilon\t\tepsilon for LJ potential\n"
                  << "\t-nu,--nu\t\tparameter for ATM potential\n"
                  << "\t-csv,--csv\t\tcsv file with particles\n"
                  << "\t-csvI,--csvInterval\t\tcsv write interval\n"
                  << "\t-blx,--blx\t\tbox length X\n"
                  << "\t-bly,--bly\t\tbox length Y\n"
                  << "\t-blZ,--blZ\t\tbox length Z\n"
                  << "\t-blcx,--blcx\t\tbottom left corner X\n"
                  << "\t-blcy,--blcy\t\tbottom left corner X\n"
                  << "\t-blcZ,--blcZ\t\tbottom left corner X\n"
                  << "\t-o,--out\t\tcsv base file name that will be used to create all csv outputs from each "
                     "simulation step\n"
                  << "-op,--outProfile\t\tjson file to save profile output\n"
                  << "-decomp,--decomp\t\t(optimal | naive) decomposition strategy for cart topology\n"
                  << "-benchyaml,--benchyaml\t\tYAML file that contains config for benchmarks\n"
                  << "\t-c,--cutoff\t\tcutoff distance" << std::endl;
    }

    cliArguments cliParse(std::vector<std::string> args) {
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
                    } else if (flag.compare("csvI") == 0 || flag.compare("-csvInterval") == 0) {
                        if (args.size() <= (size_t)(i + 1)) {
                            a.printHelp();
                            exit(1);
                        }
                        value = args[i + 1];
                        a.csvWriteInterval = std::stoi(value);
                        i++;
                    } else if (flag.compare("r") == 0 || flag.compare("-respaStepSize") == 0) {
                        if (args.size() <= (size_t)(i + 1)) {
                            a.printHelp();
                            exit(1);
                        }
                        value = args[i + 1];
                        a.respaStepSize = std::stoi(value);
                        i++;
                    } else if (flag.compare("blx") == 0 || flag.compare("-blx") == 0) {
                        if (args.size() <= (size_t)(i + 1)) {
                            a.printHelp();
                            exit(1);
                        }
                        value = args[i + 1];
                        a.boxSize[0] = std::stod(value);
                        i++;
                    } else if (flag.compare("bly") == 0 || flag.compare("-bly") == 0) {
                        if (args.size() <= (size_t)(i + 1)) {
                            a.printHelp();
                            exit(1);
                        }
                        value = args[i + 1];
                        a.boxSize[1] = std::stod(value);
                        i++;
                    } else if (flag.compare("blz") == 0 || flag.compare("-blz") == 0) {
                        if (args.size() <= (size_t)(i + 1)) {
                            a.printHelp();
                            exit(1);
                        }
                        value = args[i + 1];
                        a.boxSize[2] = std::stod(value);
                        i++;
                    } else if (flag.compare("blcx") == 0 || flag.compare("-blcx") == 0) {
                        if (args.size() <= (size_t)(i + 1)) {
                            a.printHelp();
                            exit(1);
                        }
                        value = args[i + 1];
                        a.bottomLeft[0] = std::stod(value);
                        i++;
                    } else if (flag.compare("blcy") == 0 || flag.compare("-blcy") == 0) {
                        if (args.size() <= (size_t)(i + 1)) {
                            a.printHelp();
                            exit(1);
                        }
                        value = args[i + 1];
                        a.bottomLeft[1] = std::stod(value);
                        i++;
                    } else if (flag.compare("blcz") == 0 || flag.compare("-blcz") == 0) {
                        if (args.size() <= (size_t)(i + 1)) {
                            a.printHelp();
                            exit(1);
                        }
                        value = args[i + 1];
                        a.bottomLeft[2] = std::stod(value);
                        i++;
                    } else if (flag.compare("tI") == 0 || flag.compare("-thermostatInterval") == 0) {
                        if (args.size() <= (size_t)(i + 1)) {
                            a.printHelp();
                            exit(1);
                        }
                        value = args[i + 1];
                        a.thermostatInterval = std::stoi(value);
                        a.useThermostat = true;
                        i++;
                    } else if (flag.compare("tInitT") == 0 || flag.compare("-thermostatInitialTemperature") == 0) {
                        if (args.size() <= (size_t)(i + 1)) {
                            a.printHelp();
                            exit(1);
                        }
                        value = args[i + 1];
                        a.initialTemperature = std::stod(value);
                        a.useThermostat = true;
                        i++;
                    } else if (flag.compare("tDeltaT") == 0 || flag.compare("-thermostatDeltaTemperature") == 0) {
                        if (args.size() <= (size_t)(i + 1)) {
                            a.printHelp();
                            exit(1);
                        }
                        value = args[i + 1];
                        a.deltaTemperature = std::stod(value);
                        a.useThermostat = true;
                        i++;
                    } else if (flag.compare("tTargetT") == 0 || flag.compare("-thermostatTargetTemperature") == 0) {
                        if (args.size() <= (size_t)(i + 1)) {
                            a.printHelp();
                            exit(1);
                        }
                        value = args[i + 1];
                        a.targetTemperature = std::stod(value);
                        a.useThermostat = true;
                        i++;
                    } else if (flag.compare("tAddB") == 0 || flag.compare("-thermostatAddBrownianMotion") == 0) {
                        a.addBrownianMotion = true;
                        a.useThermostat = true;
                    } else if (flag.compare("a") == 0 || flag.compare("-algorithm") == 0) {
                        if (args.size() <= (size_t)(i + 1)) {
                            a.printHelp();
                            exit(1);
                        }
                        value = args[i + 1];
                        if (value.compare("eauta") == 0) {
                            a.algorithm = AlgorithmType::EAUTAType;
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
                    } else if (flag.compare("s") == 0 || flag.compare("-sigma") == 0) {
                        if (args.size() <= (size_t)(i + 1)) {
                            a.printHelp();
                            exit(1);
                        }
                        value = args[i + 1];
                        a.sigma = std::stod(value);

                        i += 1;
                    } else if (flag.compare("e") == 0 || flag.compare("-epsilon") == 0) {
                        if (args.size() <= (size_t)(i + 1)) {
                            a.printHelp();
                            exit(1);
                        }
                        value = args[i + 1];
                        a.epsilon = std::stod(value);

                        i += 1;
                    } else if (flag.compare("nu") == 0 || flag.compare("-nu") == 0) {
                        if (args.size() <= (size_t)(i + 1)) {
                            a.printHelp();
                            exit(1);
                        }
                        value = args[i + 1];
                        a.nu = std::stod(value);

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