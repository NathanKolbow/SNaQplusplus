// Ports of SNaQ.jl's identifiability checks
// (src/network_properties/identifiability_properties.jl, `lazysnaq` branch):
// `tcgidentifiable`, `isCk`, `getdiamondns` and `getleavesbelow`.
// Used as the `--tcg_identifiable` restriction (Julia: restrictions=tcgidentifiable).
#pragma once

#include <utility>
#include <vector>

#include "network.hpp"

namespace csnaq {

// tcgidentifiable(N): true iff N's semi-directed topology and branch parameters
// are known to be identifiable (trees: true; rooted networks: false).
bool tcgidentifiable(const Network& net);

// isCk(N): (false, 0) if N is not of any class Ck, (true, k) otherwise.
std::pair<bool, int> isCk(const Network& net);

// getleavesbelow(nodes): leaves below `nodes`, following edge directions.
std::vector<NodeId> getleavesbelow(const Network& net, const std::vector<NodeId>& nodes);

// getdiamondns(cycle): the n0, n1, n2, n3 taxa counts of a 4-cycle.
std::vector<int> getdiamondns(const Network& net, const std::vector<NodeId>& cycle);

}  // namespace csnaq
