#pragma once

#include <numeric>

#include "Topology.hpp"

struct CartRank {
private:
    int dimensions;
    std::array<int, 3> rank;

public:
    CartRank() : dimensions(0), rank({0, 0, 0}) {}
    CartRank(int x) : dimensions(1), rank({x, 0, 0}) {}
    CartRank(int x, int y) : dimensions(2), rank({x, y, 0}) {}
    CartRank(int x, int y, int z) : dimensions(3), rank({x, y, z}) {}
    std::array<int, 3> GetRank() { return rank; }
    int GetDimensions() { return dimensions; }
};

class CartTopology final : public Topology {
protected:
    int dimX, dimY, dimZ;
    CartRank cartRank;
    std::vector<int> decomposition;

public:
    CartTopology(std::vector<int> decomposition);
    virtual ~CartTopology();

    CartRank GetCartRank();
    std::tuple<int, int> Shift(int dim, int dir);

    int GetLeftNeighbor(int dim);
    int GetRightNeighbor(int dim);

    void Init(std::shared_ptr<Simulation> simulation) override;

    std::array<int, 3> GetDims();
};