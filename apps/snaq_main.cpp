// `snaq`: command-line front end of SNaQ++.
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "cli.hpp"
#include "newick.hpp"
#include "readcsv.hpp"
#include "snaq.hpp"

#ifndef SNAQPP_VERSION
#define SNAQPP_VERSION "0.1.0"
#endif

int main(int argc, char** argv) {
    using namespace snaqpp;
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
            std::cout << "snaq (SNaQ++) " << SNAQPP_VERSION << "\n";
            return 0;
        }

        // Determine whether the input is a CF table (CSV) or gene trees.
        const bool use_csv = a.csv || looks_like_csv(a.input);

        if (!a.optimize.empty()) {
            std::vector<Network> nets = read_networks_arg(a.optimize);
            if (nets.size() != 1)
                throw UsageError("--optimize needs exactly one network, found " + std::to_string(nets.size()));
            if (use_csv)
                throw UsageError("--optimize with a CF table (CSV input) is not yet supported");
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
        bool single_input = start.size() == 1;

        if (use_csv) {
            // CF table input: read the CSV and pass the CFs directly to snaq.
            CsvCF csv = readtablecsv(a.input);
            const long long nq = static_cast<long long>(csv.cfs.nrow());
            std::cout << "read CF table from " << a.input << ": " << csv.taxa.size() << " taxa, " << nq
                      << " quartets" << std::endl;
            if (!csv.weights.empty()) a.opts.geneWeights = true;
            Network best = snaq_cfs(start, single_input, std::move(csv.cfs), std::move(csv.taxa), a.opts,
                                    std::move(csv.weights));
            std::cout << writenewick(best) << "\nSNaQ score = " << julia_float(best.fscore) << "\n";
        } else {
            // Gene trees input (the original path).
            std::vector<Network> genetrees = readmultinewick_file(a.input);
            if (genetrees.empty()) throw UsageError("no gene trees found in " + a.input);
            // One Newick topology behaves like snaq!(T::HybridNetwork, ...);
            // several like snaq!(Vector{HybridNetwork}, ...).
            Network best = snaq(start, single_input, genetrees, a.opts);
            std::cout << writenewick(best) << "\nSNaQ score = " << julia_float(best.fscore) << "\n";
        }
    } catch (const UsageError& e) {
        std::cerr << "ERROR: " << e.what() << "\n\n" << short_usage();
        return 1;
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
