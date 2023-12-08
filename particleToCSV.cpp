#include <getopt.h>
#include <stdio.h>

#include <array>
#include <fstream>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <tuple>

#include "ClosestPackedGenerator.hpp"
#include "ClusteredGaussGenerator.hpp"
#include "GaussGenerator.hpp"
#include "GridGenerator.hpp"
#include "ParticleGenerator.hpp"
#include "UniformGenerator.hpp"

int main(int argc, char *argv[])
{
    std::vector<std::tuple<int, double, double, double, double, double, double, double, double, double, double>>
        particles;
    int numParticles = 0;
    std::array<double, 3> distributionMean;
    std::array<double, 3> distributionStdDev;
    std::array<double, 3> bottomLeftCorner = {0.0, 0.0, 0.0};
    std::array<double, 3> boxLength;
    std::array<size_t, 3> particlesPerDim;
    int numClusters;
    double particleSpacing;
    std::string output;
    std::array<double, 3> velocity;
    double mass = 0;
    ParticleGenerator::Generator gen = ParticleGenerator::Generator::Uniform;
    uint_fast32_t seed0 = 0;
    uint_fast32_t seed1 = 0;

    static const struct option long_options[] = {{"generator", required_argument, 0, 'a'},
                                                 {"numparticles", required_argument, 0, 'b'},
                                                 {"distmean", required_argument, 0, 'c'},
                                                 {"diststddev", required_argument, 0, 'd'},
                                                 {"output", required_argument, 0, 'e'},
                                                 {"vx", required_argument, 0, 'f'},
                                                 {"vy", required_argument, 0, 'g'},
                                                 {"vz", required_argument, 0, 'h'},
                                                 {"mass", required_argument, 0, 'i'},
                                                 {"seed0", required_argument, 0, 'j'},
                                                 {"seed1", required_argument, 0, 'k'},
                                                 {"blx", required_argument, 0, 'l'},
                                                 {"bly", required_argument, 0, 'm'},
                                                 {"blz", required_argument, 0, 'n'},
                                                 {"numclusters", required_argument, 0, 'o'},
                                                 {"particlespacing", required_argument, 0, 'p'},
                                                 {"particlesperx", required_argument, 0, 'q'},
                                                 {"particlespery", required_argument, 0, 'r'},
                                                 {"particlesperz", required_argument, 0, 's'},
                                                 {0, 0, 0, 0}};

    // http://www.mario-konrad.ch/blog/programming/getopt.html
    while (1) {
        int index = -1;
        struct option *opt = 0;
        int result = getopt_long_only(argc, argv, "n:m:d:o:v:w:a:b:p:q:", long_options, &index);
        if (result == -1) break; /* end of list */
        switch (result) {
            case 'a': gen = ParticleGenerator::Str2Gen(optarg); break;
            case 'b': numParticles = std::stoi(optarg); break;
            case 'c': distributionMean[0] = distributionMean[1] = distributionMean[2] = std::stod(optarg); break;
            case 'd': distributionStdDev[0] = distributionStdDev[1] = distributionStdDev[2] = std::stod(optarg); break;
            case 'e': output = optarg; break;
            case 'f': velocity[0] = std::stod(optarg); break;
            case 'g': velocity[1] = std::stod(optarg); break;
            case 'h': velocity[2] = std::stod(optarg); break;
            case 'i': mass = std::stod(optarg); break;
            case 'j': seed0 = std::stoul(optarg); break;
            case 'k': seed1 = std::stoul(optarg); break;
            case 'l': boxLength[0] = std::stod(optarg); break;
            case 'm': boxLength[1] = std::stod(optarg); break;
            case 'n': boxLength[2] = std::stod(optarg); break;
            case 'o': numClusters = std::stoi(optarg); break;
            case 'p': particleSpacing = std::stod(optarg); break;
            case 'q': particlesPerDim[0] = std::stoul(optarg); break;
            case 'r': particlesPerDim[1] = std::stoul(optarg); break;
            case 's': particlesPerDim[2] = std::stoul(optarg); break;
            case 0:
                opt = (struct option *)&(long_options[index]);
                printf("'%s' was specified.", opt->name);
                if (opt->has_arg == required_argument) printf("Arg: <%s>", optarg);
                printf("\n");
                break;
            default: break;
        }
    }

    std::random_device rd;  // obtain a random number from hardware
    std::shared_ptr<ParticleGenerator> pGenerator(nullptr);

    switch (gen) {
        case ParticleGenerator::Generator::ClosestPacked:
            pGenerator = std::make_shared<ClosestPackedGenerator>(numParticles, velocity, boxLength, bottomLeftCorner,
                                                                  mass, seed0, seed1, particleSpacing);
            break;
        case ParticleGenerator::Generator::ClusteredGauss:
            pGenerator = std::make_shared<ClusteredGaussGenerator>(numParticles, velocity, boxLength, bottomLeftCorner,
                                                                   mass, seed0, seed1, distributionMean,
                                                                   distributionStdDev, numClusters);
            break;
        case ParticleGenerator::Generator::Gauss:
            pGenerator = std::make_shared<GaussGenerator>(numParticles, velocity, boxLength, bottomLeftCorner, mass,
                                                          seed0, seed1, distributionMean, distributionStdDev);
            break;
        case ParticleGenerator::Generator::Grid:
            pGenerator = std::make_shared<GridGenerator>(numParticles, velocity, boxLength, bottomLeftCorner, mass,
                                                         seed0, seed1, particlesPerDim, particleSpacing);
            break;
        case ParticleGenerator::Generator::Uniform:
            pGenerator = std::make_shared<UniformGenerator>(numParticles, velocity, boxLength, bottomLeftCorner, mass,
                                                            seed0, seed1);
            break;

        default: break;
    }

    pGenerator->Generate();

    particles = pGenerator->GetParticles();

    std::ofstream csvFile;
    csvFile.open(output);
    csvFile << "ID, pX, pY, pZ, vX, vY, vZ, aX, aY, aZ, m\n";
    for (size_t i = 0; i < particles.size(); i++) {
        csvFile << std::get<0>(particles[i]) << ", " << std::get<1>(particles[i]) << ", " << std::get<2>(particles[i])
                << ", " << std::get<3>(particles[i]) << ", " << std::get<4>(particles[i]) << ", "
                << std::get<5>(particles[i]) << ", " << std::get<6>(particles[i]) << ", " << std::get<7>(particles[i])
                << ", " << std::get<8>(particles[i]) << ", " << std::get<9>(particles[i]) << ", "
                << std::get<10>(particles[i]) << "\n";
    }
    csvFile.close();

    return 0;
}