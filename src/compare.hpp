// Comparing networks (PhyloNetworks compareNetworks.jl).
#pragma once

#include "network.hpp"

namespace csnaq {

// PhyloNetworks.hardwiredclusterdistance(net1, net2, rooted). With
// rooted = false and a network, every admissible root position is tried.
int hardwiredclusterdistance(const Network& net1, const Network& net2, bool rooted);

}  // namespace csnaq
