// A signature of a semidirected network's topology, used to recognize
// topologies that a search has already evaluated.
#pragma once

#include <cstdint>

#include "network.hpp"

namespace csnaq {

// Hash of the topology of `net` as a semidirected network with labeled
// leaves: tree edges are undirected, hybrid edges point to their hybrid node,
// and the position of the root, edge lengths, γ values and which hybrid edge
// is major are ignored. Internal nodes must have degree 3 (as after
// semidirectnetwork!).
//
// The hash refines node colors from the leaf labels outwards (color
// refinement). `discrete` is set to true if every node ends with its own
// color: two networks with equal signatures are then the same topology (up
// to 64-bit hash collisions). When it is false, the signature must not be
// used to conclude that two topologies are equal.
uint64_t topologysignature(const Network& net, bool* discrete);

}  // namespace csnaq
