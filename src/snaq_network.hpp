// SNaQ-specific network helpers (network_moves/misc.jl, network_utils.jl,
// network_properties/*).
#pragma once

#include "network.hpp"

namespace csnaq {

// semidirectnetwork!(N; check_conditions)
void semidirectnetwork(Network& net, bool check_conditions = true);

// isdescendantof(descendant, ancestor) (SNaQ's version, with iteration cap)
bool isdescendantof(const Network& net, NodeId descendant, NodeId ancestor);

// swapedgehybridinfo!(e1, e2)
void swapedgehybridinfo(Network& net, EdgeId e1, EdgeId e2);

// areconnected(n1, n2)
bool areconnected(const Network& net, NodeId n1, NodeId n2);

// getnetworklevel(net) and isgallednetwork(net) (SNaQ's versions)
int getnetworklevel(Network& net);
bool isgallednetwork(Network& net);

// removehybrid!(N, hyb, minor) and removeexcessreticulations!(net, hmax)
void removehybrid(Network& net, NodeId hyb, bool minor = true);
void removeexcessreticulations(Network& net, int hmax);

}  // namespace csnaq
