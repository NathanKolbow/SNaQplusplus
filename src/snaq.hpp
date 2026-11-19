// Top-level `snaq!` (inference/snaq.jl and inference/lazysnaq.jl).
#pragma once

#include <vector>

#include "network.hpp"
#include "search.hpp"

namespace snaqpp {

// Runs SNaQ from starting topologies and gene trees. `startnets` holds one
// network (single_input = true, as when `snaq!` is given a HybridNetwork) or
// one or `runs` networks (a Vector{HybridNetwork}). Writes the output files
// and returns the best network.
Network snaq(const std::vector<Network>& startnets, bool single_input, const std::vector<Network>& genetrees,
             SnaqOptions o);

// snaq! with observed CFs given directly (a non-lazy DataCF): `cfs` holds the
// CFs of all quartets of the sorted `taxa`, in SNaQ's quartet order.
// `weights` (optional): one weight per quartet, or empty (all quartets weigh 1).
Network snaq_cfs(const std::vector<Network>& startnets, bool single_input, CFMatrix cfs,
                 std::vector<std::string> taxa, SnaqOptions o, std::vector<double> weights = {});

// `snaq --optimize`: fits the branch lengths and inheritance probabilities γ
// of the fixed topology `net` to the gene trees, without changing the
// topology. Uses the quartets sampled with o.propQuartets (all quartets when
// 1, sampled with o.seed otherwise), and the parameter-optimization options of
// `o` (optimizer, opt_maxeval, ftolRel, ftolAbs, xtolRel, xtolAbs, rho,
// geneWeights, outgroup, threads). Returns the network with its score in fscore.
Network optimize_network(Network net, const std::vector<Network>& genetrees, SnaqOptions o);

}  // namespace snaqpp
