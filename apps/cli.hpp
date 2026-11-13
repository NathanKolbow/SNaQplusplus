// Command-line parsing for the `snaq` binary.
#pragma once

#include <string>
#include <vector>

#include "search.hpp"

namespace csnaq {

struct CliArgs {
    std::string input;       // -i: gene trees
    std::string start_tree;  // --start-tree: Newick string or file
    SnaqOptions opts;
    bool help = false;
    bool version = false;
};

// Throws ArgumentError on invalid usage.
CliArgs parse_cli(const std::vector<std::string>& args);
std::string usage();

}  // namespace csnaq
