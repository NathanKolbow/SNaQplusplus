// `snaq`: command-line front end of CSNaQ.
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "cli.hpp"
#include "newick.hpp"
#include "snaq.hpp"

#ifndef CSNAQ_VERSION
#define CSNAQ_VERSION "0.1.0"
#endif

int main(int argc, char** argv) {
    using namespace csnaq;
    std::vector<std::string> args(argv + 1, argv + argc);
    if (args.empty()) {
        std::cout << short_usage();
        return 0;
    }
    try {
        CliArgs a = parse_cli(args);
        if (a.help) {
            std::cout << usage();
            return 0;
        }
        if (a.version) {
            std::cout << "snaq (CSNaQ) " << CSNAQ_VERSION << "\n";
            return 0;
        }
        if (!a.optimize.empty()) {
            std::vector<Network> nets = read_networks_arg(a.optimize);
            if (nets.size() != 1)
                throw UsageError("--optimize needs exactly one network, found " + std::to_string(nets.size()));
            std::vector<Network> genetrees = readmultinewick_file(a.input);
            if (genetrees.empty()) throw UsageError("no gene trees found in " + a.input);
            Network net = optimize_network(std::move(nets[0]), genetrees, a.opts);
            const double score = net.fscore;
            const std::string nw = writenewick(net);
            if (!a.output.empty()) {
                std::ofstream f(a.output, std::ios::trunc);
                if (!f) throw std::runtime_error("cannot write " + a.output);
                f << nw << "\n";
            }
            std::cout << nw << "\nSNaQ score = " << julia_float(score) << "\n";
            return 0;
        }
        std::vector<Network> start = read_networks_arg(a.start_tree);
        if (start.empty()) throw UsageError("no starting topology found in --start-tree");
        std::vector<Network> genetrees = readmultinewick_file(a.input);
        if (genetrees.empty()) throw UsageError("no gene trees found in " + a.input);
        // One Newick topology behaves like snaq!(T::HybridNetwork, ...);
        // several like snaq!(Vector{HybridNetwork}, ...).
        bool single_input = start.size() == 1;
        Network best = snaq(start, single_input, genetrees, a.opts);
        std::cout << writenewick(best) << "\nSNaQ score = " << julia_float(best.fscore) << "\n";
    } catch (const UsageError& e) {
        std::cerr << "ERROR: " << e.what() << "\n\n" << short_usage();
        return 1;
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
