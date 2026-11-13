// Helpers shared by the C++ ports of SNaQ.jl's tests, including ports of the
// few PhyloNetworks/SNaQ functions that only the tests use
// (hardwiredclusterdistance, shrinkedge!, tcgidentifiable).
#pragma once

#include <gtest/gtest.h>

#include <initializer_list>
#include <string>
#include <vector>

#include "equations.hpp"
#include "likelihood.hpp"
#include "moves.hpp"
#include "network.hpp"
#include "network_manip.hpp"
#include "newick.hpp"
#include "snaq.hpp"
#include "snaq_network.hpp"

namespace csnaq::test {

// Path of a fixture file in tests/data.
std::string data(const std::string& file);
// Temporary directory for output files of a test.
std::string tempdir();

Network nw(const std::string& newick);
std::vector<Network> read_trees(const std::string& fixture);
std::vector<std::string> read_lines(const std::string& fixture);
CFMatrix read_cf(const std::string& fixture);

NodeId node_named(const Network& net, const std::string& name);
std::vector<NodeId> get_nodes(const Network& net, std::initializer_list<const char*> names);

// writenewick of a copy with all edge lengths missing, compared to `correct`.
bool testnewickequality(const Network& net, const std::string& correct);

// PhyloNetworks.hardwiredclusterdistance(net1, net2, rooted)
int hardwiredclusterdistance(const Network& net1, const Network& net2, bool rooted);
// PhyloNetworks.shrinkedge!(net, edge): returns true if a bad polytomy was created
bool shrinkedge(Network& net, EdgeId e);
// SNaQ.tcgidentifiable(N)
bool tcgidentifiable(const Network& N);

// A random binary tree on taxa t1..tn with internal nodes named "i<number>"
// (stand-in for SNaQ's test helper `generate_tree`).
Network random_tree(int n, Rng& rng);

// Default options for snaq! calls in tests: no output files, few threads.
SnaqOptions quiet_options();

}  // namespace csnaq::test
