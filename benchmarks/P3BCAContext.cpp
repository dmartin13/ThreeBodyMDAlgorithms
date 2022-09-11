#ifdef BENCHMARK_3BMDA

#include "P3BCAContext.hpp"

P3BCAContext::P3BCAContext(MPI_Datatype &mpiParticleType) : Context(mpiParticleType) {}

P3BCAContext::~P3BCAContext() {}

void P3BCAContext::Init(ContextArgs args)
{
    // create topology
    std::shared_ptr<CartTopology> cartTopology = std::make_shared<CartTopology>();

    // domain decomposition
    std::shared_ptr<RegularGridDecomposition> regularGridDecomposition = std::make_shared<RegularGridDecomposition>();

    // create potential
    std::shared_ptr<AxilrodTeller> axilrodTeller = std::make_shared<AxilrodTeller>(1.0);

    // create algorithm
    std::shared_ptr<P3BCA> p3bca = std::make_shared<P3BCA>(args.cutoff);

    this->simulation =
        std::make_shared<Simulation>(args.iterations, p3bca, cartTopology, axilrodTeller, regularGridDecomposition,
                                     &this->mpiParticleType, this->particles, args.deltaT, args.gForce);
}

void P3BCAContext::AfterBench(benchmark::State &state __attribute__((unused)))
{
    state.counters["cutoff_window"] =
        std::static_pointer_cast<P3BCA>(this->simulation->GetAlgorithm())->GetNumCutoffBoxes();
    state.counters["cutoff_specified"] = std::static_pointer_cast<P3BCA>(this->simulation->GetAlgorithm())->GetCutoff();
}

#endif