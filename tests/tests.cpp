#ifdef TESTS_3BMDA

#include "tests.hpp"

std::vector<Utility::Particle> particles;
MPI_Datatype mpiParticleType;
std::vector<MPI_Datatype> types;
MPI_Comm parent;

MPI_Comm interComm0;
MPI_Comm interComm1;
int myRankInterComm0;
int myRankInterComm1;

std::vector<Utility::Particle> uniformParticles;
std::vector<Utility::Particle> gridParticles;
std::vector<Utility::Particle> gaussParticles;
std::vector<Utility::Particle> closestpackedParticles;
std::vector<Utility::Particle> clusteredgaussParticles;

void generateParticles()
{
    int numParticles = 1000;
    std::array<double, 3> velocity = {0, 0, 0};
    std::array<double, 3> boxLength = {10, 10, 10};
    std::array<double, 3> bottomLeftCorner = {0, 0, 0};
    double mass = 0.0001;
    uint_fast32_t seed0 = 926762934;
    uint_fast32_t seed1 = 89347587;
    std::array<size_t, 3> particlesPerDim = {10, 10, 10};
    double particleSpacing = 0.5;
    const std::array<double, 3> distributionMean = {0.5, 0.5, 0.5};
    const std::array<double, 3> distributionStdDev = {0.5, 0.5, 0.5};
    int numClusters = 25;

    std::vector<std::tuple<int, double, double, double, double, double, double, double, double, double, double>>
        uniformParticlesTuple;
    std::vector<std::tuple<int, double, double, double, double, double, double, double, double, double, double>>
        gridParticlesTuple;
    std::vector<std::tuple<int, double, double, double, double, double, double, double, double, double, double>>
        gaussParticlesTuple;
    std::vector<std::tuple<int, double, double, double, double, double, double, double, double, double, double>>
        closestpackedParticlesTuple;
    std::vector<std::tuple<int, double, double, double, double, double, double, double, double, double, double>>
        clusteredgaussParticlesTuple;

    UniformGenerator uniformGenerator(numParticles, velocity, boxLength, bottomLeftCorner, mass, seed0, seed1);
    GridGenerator gridGenerator(numParticles, velocity, boxLength, bottomLeftCorner, mass, seed0, seed1,
                                particlesPerDim, particleSpacing);
    GaussGenerator gaussGenerator(numParticles, velocity, boxLength, bottomLeftCorner, mass, seed0, seed1,
                                  distributionMean, distributionStdDev);
    ClosestPackedGenerator closestPackedGenerator(numClusters, velocity, boxLength, bottomLeftCorner, mass, seed0,
                                                  seed1, particleSpacing);
    ClusteredGaussGenerator clusteredGaussGenerator(numParticles, velocity, boxLength, bottomLeftCorner, mass, seed0,
                                                    seed1, distributionMean, distributionStdDev, numClusters);

    uniformGenerator.Generate();
    gridGenerator.Generate();
    gaussGenerator.Generate();
    closestPackedGenerator.Generate();
    clusteredGaussGenerator.Generate();

    uniformParticlesTuple = uniformGenerator.GetParticles();
    gridParticlesTuple = gridGenerator.GetParticles();
    gaussParticlesTuple = gaussGenerator.GetParticles();
    closestpackedParticlesTuple = closestPackedGenerator.GetParticles();
    clusteredgaussParticlesTuple = clusteredGaussGenerator.GetParticles();

    Utility::getParticlesFromTuple(uniformParticlesTuple, uniformParticles);
    Utility::getParticlesFromTuple(gridParticlesTuple, gridParticles);
    Utility::getParticlesFromTuple(gaussParticlesTuple, gaussParticles);
    Utility::getParticlesFromTuple(closestpackedParticlesTuple, closestpackedParticles);
    Utility::getParticlesFromTuple(clusteredgaussParticlesTuple, clusteredgaussParticles);
}

struct CartRankTriplet {
    int a0, a1, a2;
    int b0, b1, b2;
    int c0, c1, c2;

    CartRankTriplet() : a0(0), a1(0), a2(0), b0(0), b1(0), b2(0), c0(0), c1(0), c2(0) {}
    CartRankTriplet(int a0, int a1, int a2, int b0, int b1, int b2, int c0, int c1, int c2)
        : a0(a0), a1(a1), a2(a2), b0(b0), b1(b1), b2(b2), c0(c0), c1(c1), c2(c2)
    {}

    bool operator==(const CartRankTriplet& t) const
    {
        return ((a0 == t.a0 && a1 == t.a1 && a2 == t.a2) && (b0 == t.b0 && b1 == t.b1 && b2 == t.b2) &&
                (c0 == t.c0 && c1 == t.c1 && c2 == t.c2)) ||
               ((a0 == t.a0 && a1 == t.a1 && a2 == t.a2) && (b0 == t.c0 && b1 == t.c1 && b2 == t.c2) &&
                (c0 == t.b0 && c1 == t.b1 && c2 == t.b2)) ||
               ((a0 == t.b0 && a1 == t.b1 && a2 == t.b2) && (b0 == t.a0 && b1 == t.a1 && b2 == t.a2) &&
                (c0 == t.c0 && c1 == t.c1 && c2 == t.c2)) ||
               ((a0 == t.b0 && a1 == t.b1 && a2 == t.b2) && (b0 == t.c0 && b1 == t.c1 && b2 == t.c2) &&
                (c0 == t.a0 && c1 == t.a1 && c2 == t.a2)) ||
               ((a0 == t.c0 && a1 == t.c1 && a2 == t.c2) && (b0 == t.b0 && b1 == t.b1 && b2 == t.b2) &&
                (c0 == t.a0 && c1 == t.a1 && c2 == t.a2)) ||
               ((a0 == t.c0 && a1 == t.c1 && a2 == t.c2) && (b0 == t.a0 && b1 == t.a1 && b2 == t.a2) &&
                (c0 == t.b0 && c1 == t.b1 && c2 == t.b2));
    }

    bool operator!=(const CartRankTriplet& t) const { return !(this->operator==(t)); }

    std::string toString()
    {
        return "[(" + std::to_string(a0) + ", " + std::to_string(a1) + ", " + std::to_string(a2) + "), (" +
               std::to_string(b0) + ", " + std::to_string(b1) + ", " + std::to_string(b2) + "), (" +
               std::to_string(c0) + ", " + std::to_string(c1) + ", " + std::to_string(c2) + ")]";
    }

    static MPI_Datatype GetMPIType()
    {
        // create MPI struct
        MPI_Datatype mpiTripletType;
        const int nitemsTriplet = 9;
        int blocklengthsTriplet[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
        MPI_Datatype types[9] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT};

        MPI_Aint offsetsTriplet[9];

        offsetsTriplet[0] = offsetof(CartRankTriplet, a0);
        offsetsTriplet[1] = offsetof(CartRankTriplet, a1);
        offsetsTriplet[2] = offsetof(CartRankTriplet, a2);
        offsetsTriplet[3] = offsetof(CartRankTriplet, b0);
        offsetsTriplet[4] = offsetof(CartRankTriplet, b1);
        offsetsTriplet[5] = offsetof(CartRankTriplet, b2);
        offsetsTriplet[6] = offsetof(CartRankTriplet, c0);
        offsetsTriplet[7] = offsetof(CartRankTriplet, c1);
        offsetsTriplet[8] = offsetof(CartRankTriplet, c2);

        MPI_Type_create_struct(nitemsTriplet, blocklengthsTriplet, offsetsTriplet, types, &mpiTripletType);

        return mpiTripletType;
    }
};

std::shared_ptr<Simulation> createNATAContext(int iterations, double deltaT, Eigen::Vector3d gForce)
{
    // create topology
    std::shared_ptr<RingTopology> ringTopology = std::make_shared<RingTopology>();

    // domain decomposition
    std::shared_ptr<AtomDecomposition> atomDecomposition = std::make_shared<AtomDecomposition>();

    // create potential
    std::shared_ptr<AxilrodTeller> axilrodTeller = std::make_shared<AxilrodTeller>(1.0);

    // create algorithm
    std::shared_ptr<NATA> nata = std::make_shared<NATA>();

    // set up simulation
    // int iterations, std::shared_ptr<Algorithm> algorithm, std::shared_ptr<Topology> topology,
    // std::shared_ptr<Potential> potential, std::shared_ptr<DomainDecomposition> decomposition,
    // MPI_Datatype* mpiParticleType, std::vector<Utility::Particle>& particles, double dt,
    // Eigen::Vector3d gForce
    std::shared_ptr<Simulation> simulation = std::make_shared<Simulation>(
        iterations, nata, ringTopology, axilrodTeller, atomDecomposition, &mpiParticleType, particles, deltaT, gForce);
    return simulation;
}

std::shared_ptr<Simulation> createP3BCAContext(int iterations, double deltaT, Eigen::Vector3d gForce, double cutoff,
                                               std::vector<int> decomposition)
{
    // create topology
    std::shared_ptr<CartTopology> cartTopology = std::make_shared<CartTopology>(decomposition);

    // domain decomposition
    std::shared_ptr<RegularGridDecomposition> regularGridDecomposition = std::make_shared<RegularGridDecomposition>();

    // create potential
    std::shared_ptr<AxilrodTeller> axilrodTeller = std::make_shared<AxilrodTeller>(1.0);

    // create algorithm
    std::shared_ptr<P3BCA> p3bca = std::make_shared<P3BCA>(cutoff);

    // set up simulation
    // int iterations, std::shared_ptr<Algorithm> algorithm, std::shared_ptr<Topology> topology,
    // std::shared_ptr<Potential> potential, std::shared_ptr<DomainDecomposition> decomposition,
    // MPI_Datatype* mpiParticleType, std::vector<Utility::Particle>& particles, double dt,
    // Eigen::Vector3d gForce
    std::shared_ptr<Simulation> simulation =
        std::make_shared<Simulation>(iterations, p3bca, cartTopology, axilrodTeller, regularGridDecomposition,
                                     &mpiParticleType, particles, deltaT, gForce);
    return simulation;
}

std::shared_ptr<Simulation> createAUTAContext(int iterations, double deltaT, Eigen::Vector3d gForce)
{
    // create topology
    std::shared_ptr<RingTopology> ringTopology = std::make_shared<RingTopology>();

    // domain decomposition
    std::shared_ptr<AtomDecomposition> atomDecomposition = std::make_shared<AtomDecomposition>();

    // create potential
    std::shared_ptr<AxilrodTeller> axilrodTeller = std::make_shared<AxilrodTeller>(1.0);

    // create algorithm
    std::shared_ptr<AUTA> auta = std::make_shared<AUTA>();

    // set up simulation
    // int iterations, std::shared_ptr<Algorithm> algorithm, std::shared_ptr<Topology> topology,
    // std::shared_ptr<Potential> potential, std::shared_ptr<DomainDecomposition> decomposition,
    // MPI_Datatype* mpiParticleType, std::vector<Utility::Particle>& particles, double dt,
    // Eigen::Vector3d gForce
    std::shared_ptr<Simulation> simulation = std::make_shared<Simulation>(
        iterations, auta, ringTopology, axilrodTeller, atomDecomposition, &mpiParticleType, particles, deltaT, gForce);
    return simulation;
}

int periodicDistance(int x, int y, int dim) { return std::min(abs(x - y), dim - abs(x - y)); }

Eigen::Array3i periodicDistanceA3i(Eigen::Array3i x, Eigen::Array3i y, int dim)
{
    return Eigen::Array3i(periodicDistance(x.x(), y.x(), dim), periodicDistance(x.y(), y.y(), dim),
                          periodicDistance(x.z(), y.z(), dim));
}

bool vLtS(Eigen::Array3i v, int scalar) { return (v.x() <= scalar) && (v.y() <= scalar) && (v.z() <= scalar); }

bool vLtV(Eigen::Array3i x, Eigen::Array3i y)
{
    if (x.x() != y.x()) {
        return x.x() < y.x();
    }
    if (x.y() != y.y()) {
        return x.y() < y.y();
    }
    if (x.z() != y.z()) {
        return x.z() < y.z();
    }
    return true;
}

int sgn(int value)
{
    if (value >= 0)
        return 1;
    else
        return -1;
}

int periodicDiff(int x, int y, int dim)
{
    return (abs(x - y) <= (dim / 2)) ? (x - y) : (sgn(y - x) * periodicDistance(x, y, dim));
}

Eigen::Array3i periodicDiffA3i(Eigen::Array3i x, Eigen::Array3i y, std::array<int, 3> dim)
{
    return Eigen::Array3i(periodicDiff(x.x(), y.x(), dim[0]), periodicDiff(x.y(), y.y(), dim[1]),
                          periodicDiff(x.z(), y.z(), dim[2]));
}

bool customLt(Eigen::Array3i r, Eigen::Array3i u, Eigen::Array3i v, std::array<int, 3> dim)
{
    Eigen::Vector3i diff0 = periodicDiffA3i(u, r, dim);
    Eigen::Vector3i diff1 = periodicDiffA3i(v, r, dim);
    return vLtV(diff0, diff1);
}

std::vector<Eigen::Array3i> getIntersectedCutoofWindow(std::vector<Eigen::Array3i>& a, std::vector<Eigen::Array3i>& b)
{
    std::vector<Eigen::Array3i> intersected;
    for (Eigen::Array3i& r_a : a) {
        for (Eigen::Array3i& r_b : b) {
            if ((r_a == r_b).all()) {
                intersected.push_back(r_a);
                break;
            }
        }
    }

    return intersected;
}

void p3bca(std::vector<std::tuple<std::tuple<int, int, int>, std::vector<CartRankTriplet>>>& interactions,
           std::array<int, 3> dims, std::array<int, 3> cutoffBoxes)
{
    std::vector<std::tuple<Eigen::Array3i, std::vector<Eigen::Array3i>>> U;  // all processor ranks and cutoffwindows

    // generate all processor ranks
    for (int i = 0; i < dims[0]; i++) {
        for (int j = 0; j < dims[1]; j++) {
            for (int k = 0; k < dims[2]; k++) {
                U.push_back(std::make_tuple<Eigen::Array3i, std::vector<Eigen::Array3i>>(
                    Eigen::Array3i(i, j, k), std::vector<Eigen::Array3i>()));
            }
        }
    }

    // generate all cutoff windows of all processors
    for (std::tuple<Eigen::Array3i, std::vector<Eigen::Array3i>>& r0 : U) {
        for (std::tuple<Eigen::Array3i, std::vector<Eigen::Array3i>>& r1 : U) {
            if (periodicDistance(std::get<0>(r1).x(), std::get<0>(r0).x(), dims[0]) <= cutoffBoxes[0] &&
                periodicDistance(std::get<0>(r1).y(), std::get<0>(r0).y(), dims[1]) <= cutoffBoxes[1] &&
                periodicDistance(std::get<0>(r1).z(), std::get<0>(r0).z(), dims[2]) <= cutoffBoxes[2]) {
                std::get<1>(r0).push_back(std::get<0>(r1));
            }
        }
    }

    // perform algorithm
    int rID = 0;
    int i2ID = 0;
    int i3ID = 0;
    for (std::tuple<Eigen::Array3i, std::vector<Eigen::Array3i>>& rw : U) {
        Eigen::Array3i& r = std::get<0>(rw);
        std::vector<Eigen::Array3i>& w_i = std::get<1>(rw);
        interactions.push_back(std::make_tuple(std::make_tuple(r.x(), r.y(), r.z()), std::vector<CartRankTriplet>()));
        for (Eigen::Array3i& i_2 : w_i) {
            if (customLt(r, r, i_2, dims)) {
                std::vector<Eigen::Array3i> w_2;
                for (std::tuple<Eigen::Array3i, std::vector<Eigen::Array3i>>& rw_2 : U) {
                    if ((std::get<0>(rw_2) == i_2).all()) {
                        w_2 = std::get<1>(rw_2);
                        break;
                    }
                }

                std::vector<Eigen::Array3i> w_r_i2 = getIntersectedCutoofWindow(w_i, w_2);
                for (Eigen::Array3i& i_3 : w_r_i2) {
                    if (customLt(r, i_2, i_3, dims)) {
                        // interact my particles
                        std::get<1>(interactions.back())
                            .push_back(CartRankTriplet(r.x(), r.y(), r.z(), i_2.x(), i_2.y(), i_2.z(), i_3.x(), i_3.y(),
                                                       i_3.z()));
                    }
                    i3ID++;
                }
            }
            i2ID++;
        }
        rID++;
    }
}

std::vector<Utility::Triplet> generateAllUniqueTriplets(int numProc)
{
    std::vector<Utility::Triplet> triplets;
    for (int i = 0; i < numProc; i++) {
        for (int j = i; j < numProc; j++) {
            for (int k = j; k < numProc; k++) {
                triplets.push_back(Utility::Triplet(i, j, k));
            }
        }
    }
    return triplets;
}

/**
 * @brief Test the Atom Decomposition if we have the same particles in each step without simulationstep and without
 * gravity
 *
 */
TEST(nata, test_decomposition)
{
    std::shared_ptr<Simulation> simulation = createNATAContext(0, 0.001, Eigen::Vector3d(0, 0, 0));
    simulation->Init();

    int myParticlesOldSize = simulation->GetDecomposition()->GetMyParticles().size();

    simulation->GetDecomposition()->UpdatePredictorStage(simulation->GetDeltaT());
    simulation->GetDecomposition()->Update(simulation->GetDeltaT(), simulation->GetGForce());

    int myParticlesNewSize = simulation->GetDecomposition()->GetMyParticles().size();

    GTEST_ASSERT_EQ(myParticlesOldSize, myParticlesNewSize);
}

/**
 * @brief Test the Atom Decomposition if we have the same particles in each step without simulationstep and without
 * gravity
 *
 */
TEST(auta, test_decomposition)
{
    std::shared_ptr<Simulation> simulation = createAUTAContext(0, 0.001, Eigen::Vector3d(0, 0, 0));
    simulation->Init();

    int myParticlesOldSize = simulation->GetDecomposition()->GetMyParticles().size();

    simulation->GetDecomposition()->UpdatePredictorStage(simulation->GetDeltaT());
    simulation->GetDecomposition()->Update(simulation->GetDeltaT(), simulation->GetGForce());

    int myParticlesNewSize = simulation->GetDecomposition()->GetMyParticles().size();

    GTEST_ASSERT_EQ(myParticlesOldSize, myParticlesNewSize);
}

/**
 * @brief Test the Grid Decomposition if we have the same particles in each step without simulationstep and without
 * gravity
 *
 */
TEST(p3bca, test_decomposition)
{
    int worldSize;
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

    std::shared_ptr<Simulation> simulation =
        createP3BCAContext(0, 1., Eigen::Vector3d(0, 0, 0), 0.5, Utility::getDecomposition(worldSize, decompositions));
    simulation->Init();

    int myParticlesOldSize = simulation->GetDecomposition()->GetMyParticles().size();

    simulation->GetDecomposition()->UpdatePredictorStage(simulation->GetDeltaT());
    simulation->GetDecomposition()->Update(simulation->GetDeltaT(), simulation->GetGForce());

    int myParticlesNewSize = simulation->GetDecomposition()->GetMyParticles().size();

    GTEST_ASSERT_EQ(myParticlesOldSize, myParticlesNewSize);
}

/**
 * @brief Test the Atom Decomposition if we have the same particles in each step with one simulation step
 *
 */
TEST(nata, test_decomposition_with_step)
{
    std::shared_ptr<Simulation> simulation = createNATAContext(0, 0.001, Eigen::Vector3d(0, 0, 0));
    simulation->Init();

    int myParticlesOldSize = simulation->GetDecomposition()->GetMyParticles().size();

    simulation->GetDecomposition()->UpdatePredictorStage(simulation->GetDeltaT());
    MPI_Barrier(simulation->GetTopology()->GetComm());
    simulation->GetAlgorithm()->SimulationStep();
    MPI_Barrier(simulation->GetTopology()->GetComm());
    simulation->GetDecomposition()->Update(simulation->GetDeltaT(), simulation->GetGForce());
    MPI_Barrier(simulation->GetTopology()->GetComm());

    int myParticlesNewSize = simulation->GetDecomposition()->GetMyParticles().size();

    GTEST_ASSERT_EQ(myParticlesOldSize, myParticlesNewSize);
}

/**
 * @brief Test the Atom Decomposition if we have the same particles in each step with one simulation step
 *
 */
TEST(auta, test_decomposition_with_step)
{
    std::shared_ptr<Simulation> simulation = createAUTAContext(0, 0.001, Eigen::Vector3d(0, 0, 0));
    simulation->Init();

    int myParticlesOldSize = simulation->GetDecomposition()->GetMyParticles().size();

    simulation->GetDecomposition()->UpdatePredictorStage(simulation->GetDeltaT());
    MPI_Barrier(simulation->GetTopology()->GetComm());
    simulation->GetAlgorithm()->SimulationStep();
    MPI_Barrier(simulation->GetTopology()->GetComm());
    simulation->GetDecomposition()->Update(simulation->GetDeltaT(), simulation->GetGForce());
    MPI_Barrier(simulation->GetTopology()->GetComm());

    int myParticlesNewSize = simulation->GetDecomposition()->GetMyParticles().size();

    GTEST_ASSERT_EQ(myParticlesOldSize, myParticlesNewSize);
}

/**
 * @brief Test the Grid Decomposition if we have the same particles in each step
 *
 */
TEST(p3bca, test_decomposition_with_step)
{
    int worldSize;
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

    std::shared_ptr<Simulation> simulation = createP3BCAContext(0, 1., Eigen::Vector3d(0, 0, 0), 0.5,
                                                                Utility::getDecomposition(worldSize, decompositions));
    simulation->Init();

    int myParticlesOldSize = simulation->GetDecomposition()->GetMyParticles().size();

    simulation->GetDecomposition()->UpdatePredictorStage(simulation->GetDeltaT());
    MPI_Barrier(simulation->GetTopology()->GetComm());
    simulation->GetAlgorithm()->SimulationStep();
    MPI_Barrier(simulation->GetTopology()->GetComm());
    simulation->GetDecomposition()->Update(simulation->GetDeltaT(), simulation->GetGForce());
    MPI_Barrier(simulation->GetTopology()->GetComm());

    int myParticlesNewSize = simulation->GetDecomposition()->GetMyParticles().size();

    GTEST_ASSERT_EQ(myParticlesOldSize, myParticlesNewSize);
}

TEST(nata, test_num_interactions)
{
    std::shared_ptr<Simulation> simulation = createNATAContext(0, 0.001, Eigen::Vector3d(0, 0, 0));
    simulation->Init();

    int numProcessors = simulation->GetTopology()->GetWorldSize();

    uint64_t numInteractionsExp = Utility::BinomialCoefficient(numProcessors + 2, 3) / numProcessors;

    if (numProcessors % 3 == 0) {
        // each processor calculates Utility::BinomialCoefficient(numProcessors + 2, 3) / numProcessors + 0.33 (=1)
        // interactions
        // only processeros with lowest rank calculate triplets, if they occur multiple times in one simulation step
        if (simulation->GetTopology()->GetWorldRank() / 3 == 0) {
            numInteractionsExp += 1;
        }
    }

    uint64_t numInteractionsAct = std::get<0>(simulation->GetAlgorithm()->SimulationStep());

    GTEST_ASSERT_EQ(numInteractionsExp, numInteractionsAct);
}

TEST(auta, test_num_interactions)
{
    std::shared_ptr<Simulation> simulation = createAUTAContext(0, 0.001, Eigen::Vector3d(0, 0, 0));
    simulation->Init();

    int numProcessors = simulation->GetTopology()->GetWorldSize();

    uint64_t numInteractionsExp = Utility::BinomialCoefficient(numProcessors + 2, 3) / numProcessors;

    if (numProcessors % 3 == 0) {
        // each processor calculates Utility::BinomialCoefficient(numProcessors + 2, 3) / numProcessors + 0.33 (=1)
        // interactions
        numInteractionsExp += 1;
    }

    uint64_t numInteractionsAct = std::get<0>(simulation->GetAlgorithm()->SimulationStep());

    GTEST_ASSERT_EQ(numInteractionsExp, numInteractionsAct);
}

TEST(nata, test_num_particle_interactions)
{
    std::shared_ptr<Simulation> simulation = createNATAContext(0, 0.001, Eigen::Vector3d(0, 0, 0));
    simulation->Init();

    int numParticles = simulation->GetAllParticles().size();

    uint64_t numInteractionsTotalExp = Utility::BinomialCoefficient(numParticles, 3);

    uint64_t numMyInteractionsAct = std::get<1>(simulation->GetAlgorithm()->SimulationStep());

    uint64_t numInteractionsTotalAct;

    // elements from each process are gathered in order of their rank
    MPI_Allreduce(&numMyInteractionsAct, &numInteractionsTotalAct, 1, MPI_UINT64_T, MPI_SUM,
                  simulation->GetTopology()->GetComm());

    GTEST_ASSERT_EQ(numInteractionsTotalExp, numInteractionsTotalAct);
}

TEST(auta, test_num_particle_interactions)
{
    std::shared_ptr<Simulation> simulation = createAUTAContext(0, 0.001, Eigen::Vector3d(0, 0, 0));
    simulation->Init();

    int numParticles = simulation->GetAllParticles().size();

    uint64_t numInteractionsTotalExp = Utility::BinomialCoefficient(numParticles, 3);

    uint64_t numMyInteractionsAct = std::get<1>(simulation->GetAlgorithm()->SimulationStep());

    uint64_t numInteractionsTotalAct;

    // elements from each process are gathered in order of their rank
    MPI_Allreduce(&numMyInteractionsAct, &numInteractionsTotalAct, 1, MPI_UINT64_T, MPI_SUM,
                  simulation->GetTopology()->GetComm());

    GTEST_ASSERT_EQ(numInteractionsTotalExp, numInteractionsTotalAct);
}

TEST(nata, test_processed_triplets)
{
    MPI_Datatype tripletType = Utility::Triplet::GetMPIType();
    MPI_Type_commit(&tripletType);

    std::shared_ptr<Simulation> simulation = createNATAContext(0, 0.001, Eigen::Vector3d(0, 0, 0));
    simulation->Init();

    int numProcessors = simulation->GetTopology()->GetWorldSize();

    simulation->GetAlgorithm()->SimulationStep();

    std::vector<Utility::Triplet> expectedTriplets = generateAllUniqueTriplets(numProcessors);

    std::vector<Utility::Triplet> processed = simulation->GetAlgorithm()->GetProcessed();

    std::vector<Utility::Triplet> summedProcessed;

    int numProcessed = processed.size();

    std::vector<int> allNumProcessed;

    allNumProcessed.resize(numProcessors);

    // elements from each process are gathered in order of their rank
    MPI_Allgather(&numProcessed, 1, MPI_INT, allNumProcessed.data(), 1, MPI_INT, simulation->GetTopology()->GetComm());

    std::vector<int> displacements;
    int sumDispl = 0;
    for (int i = 0; i < numProcessors; i++) {
        displacements.push_back(sumDispl);
        sumDispl += allNumProcessed[i];
    }

    summedProcessed.resize(sumDispl);

    MPI_Allgatherv(processed.data(), processed.size(), tripletType, summedProcessed.data(), allNumProcessed.data(),
                   displacements.data(), tripletType, simulation->GetTopology()->GetComm());

    bool result = true;
    for (std::vector<Utility::Triplet>::iterator it = summedProcessed.begin(); it != summedProcessed.end();) {
        std::vector<Utility::Triplet>::iterator itExp =
            std::find(expectedTriplets.begin(), expectedTriplets.end(), *it);
        if (itExp != expectedTriplets.end()) {
            it = summedProcessed.erase(it);
            expectedTriplets.erase(itExp);
        } else {
            // a triplet has not been calculated
            result = false;
            break;
        }
    }

    MPI_Type_free(&tripletType);

    // test if we have calculated all neccesary triplets, but not less
    GTEST_ASSERT_EQ(summedProcessed.size(), 0);
    GTEST_ASSERT_EQ(expectedTriplets.size(), 0);
    GTEST_ASSERT_TRUE(result);
}

TEST(auta, test_processed_triplets)
{
    MPI_Datatype tripletType = Utility::Triplet::GetMPIType();
    MPI_Type_commit(&tripletType);

    std::shared_ptr<Simulation> simulation = createAUTAContext(0, 0.001, Eigen::Vector3d(0, 0, 0));
    simulation->Init();

    int numProcessors = simulation->GetTopology()->GetWorldSize();

    simulation->GetAlgorithm()->SimulationStep();

    std::vector<Utility::Triplet> expectedTriplets = generateAllUniqueTriplets(numProcessors);

    std::vector<Utility::Triplet> processed = simulation->GetAlgorithm()->GetProcessed();

    std::vector<Utility::Triplet> summedProcessed;

    int numProcessed = processed.size();

    std::vector<int> allNumProcessed;

    allNumProcessed.resize(numProcessors);

    // elements from each process are gathered in order of their rank
    MPI_Allgather(&numProcessed, 1, MPI_INT, allNumProcessed.data(), 1, MPI_INT, simulation->GetTopology()->GetComm());

    std::vector<int> displacements;
    int sumDispl = 0;
    for (int i = 0; i < numProcessors; i++) {
        displacements.push_back(sumDispl);
        sumDispl += allNumProcessed[i];
    }

    summedProcessed.resize(sumDispl);

    MPI_Allgatherv(processed.data(), processed.size(), tripletType, summedProcessed.data(), allNumProcessed.data(),
                   displacements.data(), tripletType, simulation->GetTopology()->GetComm());

    bool result = true;
    for (std::vector<Utility::Triplet>::iterator it = summedProcessed.begin(); it != summedProcessed.end();) {
        std::vector<Utility::Triplet>::iterator itExp =
            std::find(expectedTriplets.begin(), expectedTriplets.end(), *it);
        if (itExp != expectedTriplets.end()) {
            it = summedProcessed.erase(it);
            expectedTriplets.erase(itExp);
        } else {
            // a triplet has not been calculated
            result = false;
            break;
        }
    }

    MPI_Type_free(&tripletType);

    // test if we have calculated all neccesary triplets, but not less or more
    GTEST_ASSERT_EQ(summedProcessed.size(), 0);
    GTEST_ASSERT_EQ(expectedTriplets.size(), 0);
    GTEST_ASSERT_TRUE(result);
}

TEST(p3bca, test_processed_triplets)
{
    int worldSize;
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

    bool noRedundancyInMyParticles = false;
    bool noRedundancyInAllParticles = false;
    bool correctNumberOfInteractions = false;
    bool allNecessaryTriplets = false;
    int redCount = 0;

    MPI_Datatype tripletType = Utility::Triplet::GetMPIType();
    MPI_Type_commit(&tripletType);

    std::shared_ptr<Simulation> simulation = createP3BCAContext(0, 0.001, Eigen::Vector3d(0, 0, 0), 0.5,
                                                                Utility::getDecomposition(worldSize, decompositions));

    simulation->Init();

    std::shared_ptr<RegularGridDecomposition> decomposition =
        std::static_pointer_cast<RegularGridDecomposition>(simulation->GetDecomposition());
    std::shared_ptr<P3BCA> algorithm = std::static_pointer_cast<P3BCA>(simulation->GetAlgorithm());
    std::shared_ptr<CartTopology> topology = std::static_pointer_cast<CartTopology>(simulation->GetTopology());

    simulation->GetAlgorithm()->SimulationStep();

    // calculate all triplets that have to be processed
    std::vector<std::tuple<std::tuple<int, int, int>, std::vector<CartRankTriplet>>> expectedTriplets;

    p3bca(expectedTriplets,
          std::array<int, 3>({decomposition->GetDimX(), decomposition->GetDimY(), decomposition->GetDimZ()}),
          algorithm->GetNumCutoffBoxes());

    std::vector<Utility::Triplet> processed = simulation->GetAlgorithm()->GetProcessed();

    // transform to cart coords
    std::vector<CartRankTriplet> myCartRankProcessed;
    for (Utility::Triplet& t : processed) {
        int p0Coords[3];
        int p1Coords[3];
        int p2Coords[3];
        MPI_Cart_coords(simulation->GetTopology()->GetComm(), t.a, 3, p0Coords);
        MPI_Cart_coords(simulation->GetTopology()->GetComm(), t.b, 3, p1Coords);
        MPI_Cart_coords(simulation->GetTopology()->GetComm(), t.c, 3, p2Coords);

        myCartRankProcessed.push_back(CartRankTriplet(p0Coords[0], p0Coords[1], p0Coords[2], p1Coords[0], p1Coords[1],
                                                      p1Coords[2], p2Coords[0], p2Coords[1], p2Coords[2]));
    }

    // find my subset of interactions
    std::tuple<std::tuple<int, int, int>, std::vector<CartRankTriplet>> expectedOfMyRank;
    for (auto e : expectedTriplets) {
        auto eRank = std::get<0>(e);
        std::array<int, 3UL> cR = topology->GetCartRank().GetRank();
        auto myRank = std::make_tuple(cR[0], cR[1], cR[2]);
        if (eRank == myRank) {
            expectedOfMyRank = e;
            break;
        }
    }

    // check number of interactions
    correctNumberOfInteractions = myCartRankProcessed.size() == std::get<1>(expectedOfMyRank).size();

    // check for redundancy in expectedOfMyRank
    redCount = 0;
    for (size_t i = 0; i < std::get<1>(expectedOfMyRank).size(); i++) {
        CartRankTriplet t0 = std::get<1>(expectedOfMyRank)[i];
        for (size_t j = 0; j < std::get<1>(expectedOfMyRank).size(); j++) {
            if (i != j) {
                CartRankTriplet t1 = std::get<1>(expectedOfMyRank)[j];
                if (t0 == t1) {
                    redCount++;
                }
            }
        }
    }
    noRedundancyInMyParticles = redCount == 0 ? true : false;

    // check for redundancy in all expected
    redCount = 0;
    std::vector<CartRankTriplet> allExpected;
    for (auto& e : expectedTriplets) {
        std::array<int, 3UL> cR = topology->GetCartRank().GetRank();
        auto myRank = std::make_tuple(cR[0], cR[1], cR[2]);
        if (std::get<0>(e) != myRank) {
            auto expected = std::get<1>(e);
            allExpected.insert(allExpected.end(), expected.begin(), expected.end());
        }
    }

    for (size_t i = 0; i < std::get<1>(expectedOfMyRank).size(); i++) {
        CartRankTriplet t0 = std::get<1>(expectedOfMyRank)[i];

        for (size_t j = 0; j < expectedTriplets.size(); j++) {
            auto eRank = std::get<0>(expectedTriplets[j]);
            std::array<int, 3UL> cR = topology->GetCartRank().GetRank();
            auto myRank = std::make_tuple(cR[0], cR[1], cR[2]);
            if (eRank != myRank) {
                for (auto& e : std::get<1>(expectedTriplets[j])) {
                    if (e == t0) {
                        redCount++;
                    }
                }
            }
        }
    }
    noRedundancyInAllParticles = redCount == 0 ? true : false;

    // check if we have calculated all necessary triplet and not more or less
    for (std::vector<CartRankTriplet>::iterator it0 = myCartRankProcessed.begin(); it0 != myCartRankProcessed.end();) {
        std::vector<CartRankTriplet>::iterator it1 =
            std::find(std::get<1>(expectedOfMyRank).begin(), std::get<1>(expectedOfMyRank).end(), *it0);
        if (it1 != std::get<1>(expectedOfMyRank).end()) {
            it0 = myCartRankProcessed.erase(it0);
            allNecessaryTriplets = true;
        } else {
            // a triplet has not been calculated
            allNecessaryTriplets = false;
            ++it0;
            break;
        }
    }
    if (myCartRankProcessed.size() > 0) {
        allNecessaryTriplets = false;
    }

    MPI_Type_free(&tripletType);

    GTEST_ASSERT_TRUE(allNecessaryTriplets);
    GTEST_ASSERT_TRUE(noRedundancyInMyParticles);
    GTEST_ASSERT_TRUE(noRedundancyInAllParticles);
    GTEST_ASSERT_TRUE(correctNumberOfInteractions);
}

TEST(utility, test_triplet_uniqueness)
{
    Utility::Triplet t0(1, 2, 3);
    Utility::Triplet t1(1, 3, 2);
    Utility::Triplet t2(2, 1, 3);
    Utility::Triplet t3(2, 3, 1);
    Utility::Triplet t4(3, 1, 2);
    Utility::Triplet t5(3, 2, 1);
    Utility::Triplet t6(1, 1, 3);

    GTEST_ASSERT_EQ(t0, t0);
    GTEST_ASSERT_EQ(t0, t1);
    GTEST_ASSERT_EQ(t0, t2);
    GTEST_ASSERT_EQ(t0, t3);
    GTEST_ASSERT_EQ(t0, t4);
    GTEST_ASSERT_EQ(t0, t5);
    GTEST_ASSERT_EQ(t1, t1);
    GTEST_ASSERT_EQ(t1, t2);
    GTEST_ASSERT_EQ(t1, t3);
    GTEST_ASSERT_EQ(t1, t4);
    GTEST_ASSERT_EQ(t1, t5);
    GTEST_ASSERT_EQ(t2, t2);
    GTEST_ASSERT_EQ(t2, t3);
    GTEST_ASSERT_EQ(t2, t4);
    GTEST_ASSERT_EQ(t2, t5);
    GTEST_ASSERT_EQ(t3, t3);
    GTEST_ASSERT_EQ(t3, t4);
    GTEST_ASSERT_EQ(t3, t5);
    GTEST_ASSERT_EQ(t4, t4);
    GTEST_ASSERT_EQ(t4, t5);
    GTEST_ASSERT_EQ(t5, t5);
    GTEST_ASSERT_NE(t0, t6);
    GTEST_ASSERT_NE(t1, t6);
    GTEST_ASSERT_NE(t2, t6);
    GTEST_ASSERT_NE(t3, t6);
    GTEST_ASSERT_NE(t4, t6);
    GTEST_ASSERT_NE(t5, t6);
}

TEST(utility, test_particle_constructor)
{
    Utility::Particle p0;
    Utility::Particle p1(true);
    Utility::Particle p2(0, 1., 1., 1., 2., 2., 2., 3., 3., 3., 4.);

    EXPECT_DOUBLE_EQ(p0.pX, 0.);
    EXPECT_DOUBLE_EQ(p0.pY, 0.);
    EXPECT_DOUBLE_EQ(p0.pZ, 0.);
    EXPECT_DOUBLE_EQ(p0.vX, 0.);
    EXPECT_DOUBLE_EQ(p0.vY, 0.);
    EXPECT_DOUBLE_EQ(p0.vZ, 0.);
    EXPECT_DOUBLE_EQ(p0.aX, 0.);
    EXPECT_DOUBLE_EQ(p0.aY, 0.);
    EXPECT_DOUBLE_EQ(p0.aZ, 0.);
    EXPECT_DOUBLE_EQ(p0.mass, 0.);
    EXPECT_FALSE(p0.isDummy);

    EXPECT_TRUE(p1.isDummy);

    EXPECT_DOUBLE_EQ(p2.pX, 1.);
    EXPECT_DOUBLE_EQ(p2.pY, 1.);
    EXPECT_DOUBLE_EQ(p2.pZ, 1.);
    EXPECT_DOUBLE_EQ(p2.vX, 2.);
    EXPECT_DOUBLE_EQ(p2.vY, 2.);
    EXPECT_DOUBLE_EQ(p2.vZ, 2.);
    EXPECT_DOUBLE_EQ(p2.aX, 3.);
    EXPECT_DOUBLE_EQ(p2.aY, 3.);
    EXPECT_DOUBLE_EQ(p2.aZ, 3.);
    EXPECT_DOUBLE_EQ(p2.mass, 4.);
    EXPECT_FALSE(p2.isDummy);
}

TEST(utility, test_particle_reset)
{
    Utility::Particle p;
    p.fX = p.fY = p.fZ = 0.1;
    p.ResetForce();
    EXPECT_DOUBLE_EQ(p.fX, 0.);
    EXPECT_DOUBLE_EQ(p.fY, 0.);
    EXPECT_DOUBLE_EQ(p.fZ, 0.);
}

TEST(utility, test_particle_GetR)
{
    Utility::Particle p(0, 1., 2., 3., 0., 0., 0., 0., 0., 0., 0.);
    Eigen::Vector3d r = p.GetR();
    EXPECT_DOUBLE_EQ(r.x(), 1.);
    EXPECT_DOUBLE_EQ(r.y(), 2.);
    EXPECT_DOUBLE_EQ(r.z(), 3.);
}

int main(int argc, char* argv[])
{
    int result;
    int numParticles;

    // if (argc == 1) {
    //    result = system("mpiexec -n 1 ./tests mpi");
    //} else if (strcmp(argv[1], "mpi") == 0) {

    // init MPI
    MPI_Init(&argc, &argv);

    MPI_Comm_get_parent(&parent);

    // create particleMPIType
    mpiParticleType = Utility::Particle::GetMPIType();
    MPI_Type_commit(&mpiParticleType);

    types.push_back(mpiParticleType);

    if (parent == MPI_COMM_NULL) {
        generateParticles();
        particles = uniformParticles;
        numParticles = particles.size();

        std::vector<char*> args;

        args.push_back((char*)"--gtest_color=yes");

        args.push_back((char*)"--gtest_filter=nata.*:auta.*:utility.*");

        MPI_Comm_spawn("./tests", args.data(), 16, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &interComm0, MPI_ERRCODES_IGNORE);
        MPI_Comm_rank(interComm0, &myRankInterComm0);

        // send particles to children
        MPI_Bcast(&numParticles, 1, MPI_INT, MPI_ROOT, interComm0);
        MPI_Bcast(particles.data(), numParticles, mpiParticleType, MPI_ROOT, interComm0);

        // The mpi listener executes a barrier after all test have been performed so we can run filtered tests
        // sequentially
        MPI_Barrier(interComm0);

        args.pop_back();
        args.push_back((char*)"--gtest_filter=p3bca.*");

        MPI_Comm_spawn("./tests", args.data(), 64, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &interComm1, MPI_ERRCODES_IGNORE);
        MPI_Comm_rank(interComm1, &myRankInterComm1);

        // send particles to children
        MPI_Bcast(&numParticles, 1, MPI_INT, MPI_ROOT, interComm1);
        MPI_Bcast(particles.data(), numParticles, mpiParticleType, MPI_ROOT, interComm1);

        // The mpi listener executes a barrier after all test have been performed so we can run filtered tests
        // sequentially
        MPI_Barrier(interComm1);

        MPI_Type_free(&mpiParticleType);
        MPI_Finalize();

        result = 0;
    } else {
        ::testing::InitGoogleTest(&argc, argv);

        // parse cli arguments
        std::vector<std::string> args;

        for (int i = 2; i < argc; i++) {
            args.push_back(argv[i]);
        }

        // receive particles from root
        MPI_Bcast(&numParticles, 1, MPI_INT, 0, parent);
        particles.resize(numParticles);
        MPI_Bcast(particles.data(), numParticles, mpiParticleType, 0, parent);

        // load particle input data
        // Utility::getParticlesFromCSV("tools/test3.csv", particles);

        // Add object that will finalize MPI on exit; Google Test owns this pointer
        ::testing::AddGlobalTestEnvironment(new GTestMPIListener::MPIEnvironment(parent, types));

        // Get the event listener list.
        ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();

        // Remove default listener: the default printer and the default XML printer
        ::testing::TestEventListener* l = listeners.Release(listeners.default_result_printer());

        // Adds MPI listener; Google Test owns this pointer
        listeners.Append(new GTestMPIListener::MPIWrapperPrinter(l, MPI_COMM_WORLD));

        result = RUN_ALL_TESTS();

        // finalize ... this is done by MPI listener
    }
    //} else {
    //    result = 1;
    //}

    return result;
}

#endif