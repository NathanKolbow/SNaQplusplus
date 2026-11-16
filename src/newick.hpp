// Extended Newick reading and writing (ports of PhyloNetworks' readnewick,
// readmultinewick and writenewick).
#pragma once

#include <string>
#include <vector>

#include "network.hpp"

namespace snaqpp {

// Parse one network from an extended Newick string (PhyloNetworks `readnewick`).
Network readnewick_string(const std::string& newick);

// Read every network in a file, one per line (PhyloNetworks `readmultinewick`).
// Lines that are empty or contain only whitespace are skipped.
std::vector<Network> readmultinewick_file(const std::string& path);

// Either a Newick string (if it starts with '(' after leading whitespace) or
// the path of a file with one or more Newick lines.
std::vector<Network> read_networks_arg(const std::string& arg);

// PhyloNetworks `writenewick(net; round, digits, di, internallabel)`.
// May modify `net` (it calls directedges!, as PhyloNetworks does).
std::string writenewick(Network& net, bool round = false, int digits = 3, bool di = false,
                        bool internallabel = true);

}  // namespace snaqpp
