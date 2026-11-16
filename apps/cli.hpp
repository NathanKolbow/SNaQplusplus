// Command-line parsing for the `snaq` binary.
#pragma once

#include <string>
#include <vector>

#include "search.hpp"

namespace snaqpp {

struct CliArgs {
    std::string input;       // -i: gene trees
    std::string start_tree;  // --start-tree: Newick string or file
    std::string optimize;    // --optimize: network file whose parameters alone are optimized
    std::string output;      // with --optimize: -o, the file the optimized network is written to ("": none)
    SnaqOptions opts;
    bool help = false;
    bool version = false;
};

// Invalid command-line usage (printed together with the USAGE message).
struct UsageError : public ArgumentError {
    explicit UsageError(const std::string& msg) : ArgumentError(msg) {}
};

// Throws UsageError on invalid usage.
CliArgs parse_cli(const std::vector<std::string>& args);
std::string usage();
// Short USAGE message with the most common flags.
std::string short_usage();

}  // namespace snaqpp
