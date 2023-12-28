#include "CartTopology.hpp"

CartTopology::CartTopology(std::vector<int> decomposition) : decomposition(decomposition) {}
CartTopology::~CartTopology() {}

void CartTopology::Init(std::shared_ptr<Simulation> simulation) {
    Topology::Init(simulation);

    int numDims = this->decomposition.size();

    this->dimX = this->decomposition[0];
    this->dimY = 1;
    this->dimZ = 1;
    if (numDims > 1) {
        dimY = this->decomposition[1];
    }
    if (numDims > 2) {
        dimZ = this->decomposition[2];
    }

    if (numDims == 1) {
        int dims[1] = {dimX};
        int periods[1] = {1};

        MPI_Cart_create(MPI_COMM_WORLD, 1, dims, periods, 1, &this->comm);
        MPI_Comm_size(this->comm, &this->worldSize);
        MPI_Comm_rank(this->comm, &this->worldRank);

        int coords[1];
        MPI_Cart_coords(this->comm, this->worldRank, 1, coords);
        this->cartRank = CartRank(coords[0]);
    } else if (numDims == 2) {
        int dims[2] = {dimX, dimY};
        int periods[2] = {1, 1};

        MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &this->comm);
        MPI_Comm_size(this->comm, &this->worldSize);
        MPI_Comm_rank(this->comm, &this->worldRank);

        int coords[2];
        MPI_Cart_coords(this->comm, this->worldRank, 2, coords);
        this->cartRank = CartRank(coords[0], coords[1]);
    } else if (numDims == 3) {
        int dims[3] = {dimX, dimY, dimZ};
        int periods[3] = {1, 1, 1};

        MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 1, &this->comm);
        MPI_Comm_size(this->comm, &this->worldSize);
        MPI_Comm_rank(this->comm, &this->worldRank);

        int coords[3];
        MPI_Cart_coords(this->comm, this->worldRank, 3, coords);
        this->cartRank = CartRank(coords[0], coords[1], coords[2]);
    } else {
        exit(1);
    }
}

std::tuple<int, int> CartTopology::Shift(int dim, int dir) {
    int source, dest;
    MPI_Cart_shift(this->comm, dim, dir, &source, &dest);
    return std::tuple<int, int>{source, dest};
}

int CartTopology::GetLeftNeighbor(int dim) { return std::get<0>(this->Shift(dim, 1)); }
int CartTopology::GetRightNeighbor(int dim) { return std::get<1>(this->Shift(dim, 1)); }

CartRank CartTopology::GetCartRank() { return this->cartRank; }

std::array<int, 3> CartTopology::GetDims() { return std::array<int, 3>({this->dimX, this->dimY, this->dimZ}); }