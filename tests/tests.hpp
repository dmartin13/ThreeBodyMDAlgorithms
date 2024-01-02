#pragma once

#ifdef TESTS_3BMDA

#include <gtest/gtest.h>
#include <mpi.h>

#include <Eigen/Dense>

#include "../tools/ClosestPackedGenerator.hpp"
#include "../tools/ClusteredGaussGenerator.hpp"
#include "../tools/GaussGenerator.hpp"
#include "../tools/GridGenerator.hpp"
#include "../tools/ParticleGenerator.hpp"
#include "../tools/UniformGenerator.hpp"
#include "algorithm/Algorithm.hpp"
#include "decomposition/AtomDecomposition.hpp"
#include "decomposition/DomainDecomposition.hpp"
#include "fwd.hpp"
#include "gtest_mpi_listener.hpp"
#include "helpers/axilrodTellerMuto.hpp"
#include "helpers/lennardJones.hpp"
#include "potential/AxilrodTeller.hpp"
#include "potential/LennardJones.hpp"
#include "potential/Potential.hpp"
#include "simulation/Simulation.hpp"
#include "topology/RingTopology.hpp"
#include "topology/Topology.hpp"
#include "utility/utility.hpp"

#endif