// Top-level `snaq!` (inference/snaq.jl and inference/lazysnaq.jl).
#pragma once

#include <vector>

#include "network.hpp"
#include "search.hpp"

namespace csnaq {

// Runs SNaQ from starting topologies and gene trees. `startnets` holds one
// network (single_input = true, as when `snaq!` is given a HybridNetwork) or
// one or `runs` networks (a Vector{HybridNetwork}). Writes the output files
// and returns the best network.
Network snaq(const std::vector<Network>& startnets, bool single_input, const std::vector<Network>& genetrees,
             SnaqOptions o);

// snaq! with observed CFs given directly (a non-lazy DataCF): `cfs` holds the
// CFs of all quartets of the sorted `taxa`, in SNaQ's quartet order.
Network snaq_cfs(const std::vector<Network>& startnets, bool single_input, CFMatrix cfs,
                 std::vector<std::string> taxa, SnaqOptions o);

}  // namespace csnaq
