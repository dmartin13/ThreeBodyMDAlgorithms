#include "RingTopology.hpp"

RingTopology::RingTopology() {}
RingTopology::~RingTopology() {}

int RingTopology::GetLeftNeighbor() { return this->leftNeighbor; }
int RingTopology::GetRightNeighbor() { return this->rightNeighbor; }

void RingTopology::Init(std::shared_ptr<Simulation> simulation) {
    Topology::Init(simulation);

    MPI_Comm_size(MPI_COMM_WORLD, &this->worldSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &this->worldRank);

    // establish ring communication structure
    this->leftNeighbor = Utility::mod((this->worldRank - 1), this->worldSize);
    this->rightNeighbor = (this->worldRank + 1) % this->worldSize;

    this->comm = MPI_COMM_WORLD;
}