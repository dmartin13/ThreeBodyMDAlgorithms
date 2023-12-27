#include "AxilrodTeller.hpp"

AxilrodTeller::AxilrodTeller(double nu) : nu(nu) {}

void AxilrodTeller::CalculateForces(Utility::Particle &i, Utility::Particle &j, Utility::Particle &k)
{
    /*
    #ifdef PROFILE_3BMDA
        std::chrono::time_point<std::chrono::system_clock> start;
        std::chrono::time_point<std::chrono::system_clock> end;
        start = std::chrono::system_clock::now();
    #endif
    */

    // see "The Role of Three-Body Interactions on the Equilibrium and Non-Equilibrium Properties of Fluids from
    // Molecular Simulation" for that

    // we assume forces are set to 0 for each particle

    // sides of triangle
    Eigen::Vector3d a = (j.GetR() - i.GetR());
    Eigen::Vector3d b = (i.GetR() - k.GetR());
    Eigen::Vector3d c = (j.GetR() - k.GetR());

    // distances and powers of distances

    double d1a = a.norm();
    double d1b = b.norm();
    double d1c = c.norm();

    double d2a = d1a * d1a;
    double d2b = d1b * d1b;
    double d2c = d1c * d1c;

    double d3a = d2a * d1a;
    double d3b = d2b * d1b;
    double d3c = d2c * d1c;

    double d4a = d3a * d1a;
    double d4b = d3b * d1b;
    double d4c = d3c * d1c;

    double d5a = d4a * d1a;
    double d5b = d4b * d1b;
    double d5c = d4c * d1c;

    double d6a = d5a * d1a;
    double d6b = d5b * d1b;
    double d6c = d5c * d1c;

    double dXa = j.pX - i.pX;
    double dYa = j.pY - i.pY;
    double dZa = j.pZ - i.pZ;

    double dXb = k.pX - i.pX;
    double dYb = k.pY - i.pY;
    double dZb = k.pZ - i.pZ;

    double dXc = k.pX - j.pX;
    double dYc = k.pY - j.pY;
    double dZc = k.pZ - j.pZ;

    double dVdRa, dVdRb, dVdRc;

    // this is the gradient of the axilrodteller potential
    dVdRa = (3. * this->nu / (8. * d1a)) *
            (-8. / (d4a * d3b * d3c) - 1. / (d5b * d5c) + 5. * d1b / (d6a * d5c) + 5. * d1c / (d6a * d5b) -
             1. / (d2a * d3b * d5c) - 1. / (d2a * d5b * d3c) - 3. / (d4a * d1b * d5c) - 3. / (d4a * d5b * d1c) -
             5. / (d6a * d1b * d3c) - 5. / (d6a * d3b * d1c) + 6. / (d4a * d3b * d3c));

    dVdRb = (3. * this->nu / (8. * d1b)) *
            (-8. / (d4b * d3a * d3c) - 1. / (d5a * d5c) + 5. * d1a / (d6b * d5c) + 5. * d1c / (d6b * d5a) -
             1. / (d2b * d3a * d5c) - 1. / (d2b * d5a * d3c) - 3. / (d4b * d1a * d5c) - 3. / (d4b * d5a * d1c) -
             5. / (d6b * d1a * d3c) - 5. / (d6b * d3a * d1c) + 6. / (d4b * d3a * d3c));

    dVdRc = (3. * this->nu / (8. * d1c)) *
            (-8. / (d4c * d3b * d3a) - 1. / (d5b * d5a) + 5. * d1b / (d6c * d5a) + 5. * d1a / (d6c * d5b) -
             1. / (d2c * d3b * d5a) - 1. / (d2c * d5b * d3a) - 3. / (d4c * d1b * d5a) - 3. / (d4c * d5b * d1a) -
             5. / (d6c * d1b * d3a) - 5. / (d6c * d3b * d1a) + 6. / (d4c * d3b * d3a));

    double newIfXContrib = dXa * dVdRa + dXb * dVdRb;
    double newIfYContrib = dYa * dVdRa + dYb * dVdRb;
    double newIfZContrib = dZa * dVdRa + dZb * dVdRb;

    double newJfXContrib = dXa * (-dVdRa) + dXc * dVdRc;
    double newJfYContrib = dYa * (-dVdRa) + dYc * dVdRc;
    double newJfZContrib = dZa * (-dVdRa) + dZc * dVdRc;

    double newKfXContrib = dXb * (-dVdRb) + dXc * (-dVdRc);
    double newKfYContrib = dYb * (-dVdRb) + dYc * (-dVdRc);
    double newKfZContrib = dZb * (-dVdRb) + dZc * (-dVdRc);

#if defined(USE_OMP) && defined(OPENMPAVAIL)
    {
#pragma omp atomic
        i.fX -= newIfXContrib;
#pragma omp atomic
        i.fY -= newIfYContrib;
#pragma omp atomic
        i.fZ -= newIfZContrib;

#pragma omp atomic
        j.fX -= newJfXContrib;
#pragma omp atomic
        j.fY -= newJfYContrib;
#pragma omp atomic
        j.fZ -= newJfZContrib;

#pragma omp atomic
        k.fX -= newKfXContrib;
#pragma omp atomic
        k.fY -= newKfYContrib;
#pragma omp atomic
        k.fZ -= newKfZContrib;
    }
#else
    i.f1X -= newIfXContrib;
    i.f1Y -= newIfYContrib;
    i.f1Z -= newIfZContrib;

    j.f1X -= newJfXContrib;
    j.f1Y -= newJfYContrib;
    j.f1Z -= newJfZContrib;

    k.f1X -= newKfXContrib;
    k.f1Y -= newKfYContrib;
    k.f1Z -= newKfZContrib;
#endif
    /*
    #ifdef PROFILE_3BMDA
        end = std::chrono::system_clock::now();
        auto elapsed_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

        // check for over & underflow. https://stackoverflow.com/a/1514309
        if (elapsed_time.count() > 0 && this->timeAcc > std::numeric_limits<int64_t>::max() - elapsed_time.count()) {
            std::cout << "Overflow Warning for profiling in CalculateForces" << std::endl;
        }
        if (elapsed_time.count() < 0 && this->timeAcc < std::numeric_limits<int64_t>::min() - elapsed_time.count()) {
            std::cout << "Underflow Warning for profiling in CalculateForces" << std::endl;
        }

        this->timeAcc += elapsed_time.count();

        this->counter++;
    #endif
    */
}