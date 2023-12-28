#pragma once

#include <mpi.h>

#include <Eigen/Dense>
#include <iostream>
#include <string>

namespace Utility {
    struct Particle {
        int ID;
        double pX, pY, pZ;
        double vX, vY, vZ;
        double oldF0X, oldF0Y, oldF0Z;
        double f0X, f0Y, f0Z;
        double f1X, f1Y, f1Z;
        double mass;
        bool isDummy;

        Particle()
            : ID(-1), pX(0.0), pY(0.0), pZ(0.0), vX(0.0), vY(0.0), vZ(0.0), oldF0X(0.0), oldF0Y(0.0), oldF0Z(0.0),
              f0X(0.0), f0Y(0.0), f0Z(0.0), f1X(0.0), f1Y(0.0), f1Z(0.0), mass(0), isDummy(false) {}
        Particle(bool isDummy) : isDummy(isDummy) {}
        Particle(int ID, double pX, double pY, double pZ, double vX, double vY, double vZ, double mass)
            : ID(ID), pX(pX), pY(pY), pZ(pZ), vX(vX), vY(vY), vZ(vZ), oldF0X(0.0), oldF0Y(0.0), oldF0Z(0.0), f0X(0.0),
              f0Y(0.0), f0Z(0.0), f1X(0.0), f1Y(0.0), f1Z(0.0), mass(mass), isDummy(false) {}

        std::string toString() {
            std::string result = std::to_string(ID) + ": ";
            result.append("(" + std::to_string(pX) + ", " + std::to_string(pY) + ", " + std::to_string(pZ) +
                          ", isDummy: " + (isDummy ? "true" : "false") + ")");
            return result;
        }

        double GetSqrDist(Particle o) {
            return (o.pX - pX) * (o.pX - pX) + (o.pY - pY) * (o.pY - pY) + (o.pZ - pZ) * (o.pZ - pZ);
        }

        double GetSqrDistPeriodic(Particle o, Eigen::Array3d physicalDomainSize) {
            double absXDist = std::abs(o.pX - pX);
            double absYDist = std::abs(o.pY - pY);
            double absZDist = std::abs(o.pZ - pZ);
            double xDistP = std::min(absXDist, physicalDomainSize[0] - absXDist);
            double yDistP = std::min(absYDist, physicalDomainSize[1] - absYDist);
            double zDistP = std::min(absZDist, physicalDomainSize[2] - absZDist);
            return xDistP * xDistP + yDistP * yDistP + zDistP * zDistP;
        }

        double GetDist(Particle o) {
            return std::sqrt((o.pX - pX) * (o.pX - pX) + (o.pY - pY) * (o.pY - pY) + (o.pZ - pZ) * (o.pZ - pZ));
        }

        Eigen::Array3d GetR() { return Eigen::Array3d(pX, pY, pZ); }

        static MPI_Datatype GetMPIType() {
            // create MPI struct
            MPI_Datatype mpiParticleType;
            const int nitemsParticle = 18;
            int blocklengthsParticle[18] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
            MPI_Datatype types[18] = {MPI_INT,    MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
                                      MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
                                      MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_C_BOOL};

            MPI_Aint offsetsParticle[18];

            offsetsParticle[0] = offsetof(Utility::Particle, ID);
            offsetsParticle[1] = offsetof(Utility::Particle, pX);
            offsetsParticle[2] = offsetof(Utility::Particle, pY);
            offsetsParticle[3] = offsetof(Utility::Particle, pZ);
            offsetsParticle[4] = offsetof(Utility::Particle, vX);
            offsetsParticle[5] = offsetof(Utility::Particle, vY);
            offsetsParticle[6] = offsetof(Utility::Particle, vZ);
            offsetsParticle[7] = offsetof(Utility::Particle, oldF0X);
            offsetsParticle[8] = offsetof(Utility::Particle, oldF0Y);
            offsetsParticle[9] = offsetof(Utility::Particle, oldF0Z);
            offsetsParticle[10] = offsetof(Utility::Particle, f0X);
            offsetsParticle[11] = offsetof(Utility::Particle, f0Y);
            offsetsParticle[12] = offsetof(Utility::Particle, f0Z);
            offsetsParticle[13] = offsetof(Utility::Particle, f1X);
            offsetsParticle[14] = offsetof(Utility::Particle, f1Y);
            offsetsParticle[15] = offsetof(Utility::Particle, f1Z);
            offsetsParticle[16] = offsetof(Utility::Particle, mass);
            offsetsParticle[17] = offsetof(Utility::Particle, isDummy);

            MPI_Type_create_struct(nitemsParticle, blocklengthsParticle, offsetsParticle, types, &mpiParticleType);

            return mpiParticleType;
        }
    };

    struct Triplet {
        int a, b, c;

        Triplet() : a(0), b(0), c(0) {}
        Triplet(int a, int b, int c) : a(a), b(b), c(c) {}

        bool operator==(const Triplet& t) const {
            return (a == t.a && b == t.b && c == t.c) || (a == t.a && c == t.b && b == t.c) ||
                   (b == t.a && a == t.b && c == t.c) || (b == t.a && c == t.b && a == t.c) ||
                   (c == t.a && a == t.b && b == t.c) || (c == t.a && b == t.b && a == t.c);
        }

        bool operator!=(const Triplet& t) const { return !(this->operator==(t)); }

        std::string toString() {
            return "(" + std::to_string(a) + ", " + std::to_string(b) + ", " + std::to_string(c) + ")";
        }

        static MPI_Datatype GetMPIType() {
            // create MPI struct
            MPI_Datatype mpiTripletType;
            const int nitemsTriplet = 3;
            int blocklengthsTriplet[3] = {1, 1, 1};
            MPI_Datatype types[3] = {MPI_INT, MPI_INT, MPI_INT};

            MPI_Aint offsetsTriplet[3];

            offsetsTriplet[0] = offsetof(Utility::Triplet, a);
            offsetsTriplet[1] = offsetof(Utility::Triplet, b);
            offsetsTriplet[2] = offsetof(Utility::Triplet, c);

            MPI_Type_create_struct(nitemsTriplet, blocklengthsTriplet, offsetsTriplet, types, &mpiTripletType);

            return mpiTripletType;
        }
    };
}  // namespace Utility
