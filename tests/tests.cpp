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

void generateParticles() {
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

    std::vector<std::tuple<int, double, double, double, double, double, double, double>> uniformParticlesTuple;
    std::vector<std::tuple<int, double, double, double, double, double, double, double>> gridParticlesTuple;
    std::vector<std::tuple<int, double, double, double, double, double, double, double>> gaussParticlesTuple;
    std::vector<std::tuple<int, double, double, double, double, double, double, double>> closestpackedParticlesTuple;
    std::vector<std::tuple<int, double, double, double, double, double, double, double>> clusteredgaussParticlesTuple;

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
        : a0(a0), a1(a1), a2(a2), b0(b0), b1(b1), b2(b2), c0(c0), c1(c1), c2(c2) {}

    bool operator==(const CartRankTriplet& t) const {
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

    std::string toString() {
        return "[(" + std::to_string(a0) + ", " + std::to_string(a1) + ", " + std::to_string(a2) + "), (" +
               std::to_string(b0) + ", " + std::to_string(b1) + ", " + std::to_string(b2) + "), (" +
               std::to_string(c0) + ", " + std::to_string(c1) + ", " + std::to_string(c2) + ")]";
    }

    static MPI_Datatype GetMPIType() {
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

int periodicDistance(int x, int y, int dim) { return std::min(abs(x - y), dim - abs(x - y)); }

Eigen::Array3i periodicDistanceA3i(Eigen::Array3i x, Eigen::Array3i y, int dim) {
    return Eigen::Array3i(periodicDistance(x.x(), y.x(), dim), periodicDistance(x.y(), y.y(), dim),
                          periodicDistance(x.z(), y.z(), dim));
}

bool vLtS(Eigen::Array3i v, int scalar) { return (v.x() <= scalar) && (v.y() <= scalar) && (v.z() <= scalar); }

bool vLtV(Eigen::Array3i x, Eigen::Array3i y) {
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

int sgn(int value) {
    if (value >= 0)
        return 1;
    else
        return -1;
}

int periodicDiff(int x, int y, int dim) {
    return (abs(x - y) <= (dim / 2)) ? (x - y) : (sgn(y - x) * periodicDistance(x, y, dim));
}

Eigen::Array3i periodicDiffA3i(Eigen::Array3i x, Eigen::Array3i y, std::array<int, 3> dim) {
    return Eigen::Array3i(periodicDiff(x.x(), y.x(), dim[0]), periodicDiff(x.y(), y.y(), dim[1]),
                          periodicDiff(x.z(), y.z(), dim[2]));
}

bool customLt(Eigen::Array3i r, Eigen::Array3i u, Eigen::Array3i v, std::array<int, 3> dim) {
    Eigen::Vector3i diff0 = periodicDiffA3i(u, r, dim);
    Eigen::Vector3i diff1 = periodicDiffA3i(v, r, dim);
    return vLtV(diff0, diff1);
}

std::vector<Eigen::Array3i> getIntersectedCutoofWindow(std::vector<Eigen::Array3i>& a, std::vector<Eigen::Array3i>& b) {
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

std::vector<Utility::Triplet> generateAllUniqueTriplets(int numProc) {
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

TEST(utility, test_triplet_uniqueness) {
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

TEST(utility, test_particle_constructor) {
    Utility::Particle p0;
    Utility::Particle p1(true);
    Utility::Particle p2(0, 1., 1., 1., 2., 2., 2., 4.);

    EXPECT_DOUBLE_EQ(p0.pX, 0.);
    EXPECT_DOUBLE_EQ(p0.pY, 0.);
    EXPECT_DOUBLE_EQ(p0.pZ, 0.);
    EXPECT_DOUBLE_EQ(p0.vX, 0.);
    EXPECT_DOUBLE_EQ(p0.vY, 0.);
    EXPECT_DOUBLE_EQ(p0.vZ, 0.);
    EXPECT_DOUBLE_EQ(p0.mass, 0.);
    EXPECT_FALSE(p0.isDummy);

    EXPECT_TRUE(p1.isDummy);

    EXPECT_DOUBLE_EQ(p2.pX, 1.);
    EXPECT_DOUBLE_EQ(p2.pY, 1.);
    EXPECT_DOUBLE_EQ(p2.pZ, 1.);
    EXPECT_DOUBLE_EQ(p2.vX, 2.);
    EXPECT_DOUBLE_EQ(p2.vY, 2.);
    EXPECT_DOUBLE_EQ(p2.vZ, 2.);
    EXPECT_DOUBLE_EQ(p2.mass, 4.);
    EXPECT_FALSE(p2.isDummy);
}

TEST(utility, test_particle_GetR) {
    Utility::Particle p(0, 1., 2., 3., 0., 0., 0., 0.);
    Eigen::Vector3d r = p.GetR();
    EXPECT_DOUBLE_EQ(r.x(), 1.);
    EXPECT_DOUBLE_EQ(r.y(), 2.);
    EXPECT_DOUBLE_EQ(r.z(), 3.);
}

int main(int argc, char* argv[]) {
    int result;
    int numParticles;

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

        // args.push_back((char*)"--gtest_filter=nata.*:auta.*:utility.*");

        MPI_Comm_spawn("./tests", args.data(), 16, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &interComm0, MPI_ERRCODES_IGNORE);
        MPI_Comm_rank(interComm0, &myRankInterComm0);

        // send particles to children
        MPI_Bcast(&numParticles, 1, MPI_INT, MPI_ROOT, interComm0);
        MPI_Bcast(particles.data(), numParticles, mpiParticleType, MPI_ROOT, interComm0);

        // The mpi listener executes a barrier after all test have been performed so we can run filtered tests
        // sequentially
        MPI_Barrier(interComm0);

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

    return result;
}

#endif