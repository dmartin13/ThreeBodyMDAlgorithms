#pragma once

#ifdef BENCHMARK_3BMDA

#include "Context.hpp"

class AUTAContext : public Context {
private:
public:
    AUTAContext(MPI_Datatype &mpiParticleType);
    ~AUTAContext();

    void Init(ContextArgs args) override;
    void AfterBench(benchmark::State &state) override;
};

#endif