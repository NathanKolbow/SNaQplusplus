// Benchmark of numerical parameter optimization (fitnumericalparameters!):
// fits the branch lengths and γ of a fixed network to the observed CFs of a
// set of gene trees and reports the SNaQ score and the time taken. Compared
// against SNaQ.jl by tools/julia/bench_fit.jl.
//
// usage: csnaq_bench_fit <network newick file> <gene trees file> <maxeval> <threads> [repeats]
// prints: ntaxa nhybrids nquartets nparams score fit_seconds equations_seconds peak_rss_MB
#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>

#include "equations.hpp"
#include "newick.hpp"
#include "optimize.hpp"
#include "quartets.hpp"
#include "snaq_network.hpp"

using namespace csnaq;

int main(int argc, char** argv) {
    if (argc < 5) {
        std::cerr << "usage: csnaq_bench_fit <network> <gene trees> <maxeval> <threads> [repeats]\n";
        return 2;
    }
    try {
        Network net0 = read_networks_arg(argv[1]).at(0);
        std::vector<Network> trees = readmultinewick_file(argv[2]);
        int maxeval = std::stoi(argv[3]);
        int nthreads = std::stoi(argv[4]);
        int repeats = argc > 5 ? std::stoi(argv[5]) : 3;

        ObservedCFs obs(trees, sorted_union_taxa(trees), false);
        std::vector<long long> ranks(static_cast<size_t>(obs.nquartets()));
        std::iota(ranks.begin(), ranks.end(), 0LL);
        CFMatrix q;
        q.v.resize(3 * ranks.size());
        obs.compute_rows(ranks, q.v.data(), nthreads);

        double best_fit = 1e300, best_eq = 1e300, score = 0.0;
        size_t np = 0;
        Network net;
        for (int r = 0; r < repeats; ++r) {
            net = net0.deepcopynetwork();
            semidirectnetwork(net);
            for (EdgeId e : net.edge)
                if (!net.N(net.getchild(e)).leaf) net.E(e).length = 1.0;  // common starting values
            for (NodeId h : net.hybrid) {
                net.E(net.getparentedge(h)).gamma = 0.7;
                net.E(net.getparentedgeminor(h)).gamma = 0.3;
            }
            double t0 = wall_time();
            EquationSet es = findquartetequations(net, nthreads);
            double t1 = wall_time();
            FitOptions fo;
            fo.maxeval = maxeval;
            score = fitnumericalparameters(net, es.eqns, q, 0.0, fo, nthreads);
            double t2 = wall_time();
            best_eq = std::min(best_eq, t1 - t0);
            best_fit = std::min(best_fit, t2 - t1);
            np = es.info.params.size();
        }
        std::cout.precision(17);
        std::cout << net.numtaxa << " " << net.numhybrids << " " << obs.nquartets() << " " << np << " " << score << " "
                  << best_fit << " " << best_eq << " " << static_cast<double>(peak_rss_bytes()) / (1 << 20) << "\n";
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
