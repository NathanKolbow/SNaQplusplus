// Accuracy evaluation for the benchmark framework (benchmarks/bench.py).
//
// usage: csnaq_eval <estimated network> <true network> <gene trees file> [threads]
//   <estimated network> and <true network>: a Newick string or a file whose
//   first line is a Newick network; for a SNaQ ".out" file, the first line is
//   used up to its first ';'.
//
// Prints one JSON object with:
//   hwcd            unrooted hardwired-cluster distance between the estimated
//                   and the true network (PhyloNetworks' hardwiredclusterdistance)
//   score_estimated SNaQ score of the estimated network (with its parameters)
//                   on the observed CFs of all quartets of the gene trees
//   score_truth     same for the true network (fitted on the data: its
//                   parameters are optimized first, so the two are comparable)
//   nhybrids        number of reticulations of the estimated network
#include <cstdio>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>

#include "compare.hpp"
#include "likelihood.hpp"
#include "newick.hpp"
#include "optimize.hpp"
#include "quartets.hpp"
#include "snaq_network.hpp"

using namespace csnaq;

static Network read_first(const std::string& arg) {
    std::string s = arg;
    size_t p = s.find_first_not_of(" \t\r\n");
    if (p == std::string::npos || s[p] != '(') {
        std::ifstream in(arg);
        if (!in) throw SnaqError("cannot open " + arg);
        std::getline(in, s);
    }
    size_t semi = s.find(';');
    if (semi == std::string::npos) throw SnaqError("no Newick network in " + arg);
    return readnewick_string(s.substr(0, semi + 1));
}

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cerr << "usage: csnaq_eval <estimated network> <true network> <gene trees file> [threads]\n";
        return 2;
    }
    try {
        Network est = read_first(argv[1]);
        Network truth = read_first(argv[2]);
        std::vector<Network> gts = readmultinewick_file(argv[3]);
        int nthreads = argc > 4 ? std::stoi(argv[4]) : 1;

        int hwcd = hardwiredclusterdistance(est, truth, false);

        ObservedCFs obs(gts, sorted_union_taxa(gts), false);
        std::vector<long long> ranks(static_cast<size_t>(obs.nquartets()));
        std::iota(ranks.begin(), ranks.end(), 0LL);
        CFMatrix q;
        q.v.resize(3 * ranks.size());
        obs.compute_rows(ranks, q.v.data(), nthreads);

        // leaf edges may be missing (-1) in output files; they do not matter
        double s_est = computeSNaQscore_net(est, q, 0.0, nthreads);
        Network t = truth.deepcopynetwork();
        semidirectnetwork(t);
        FitOptions fo;
        fo.maxeval = 1000;
        fo.method = FitMethod::newton;
        double s_truth = fitnumericalparameters_all(t, q, 0.0, fo, nthreads);

        // the estimated topology with its parameters fitted to all quartets, like the true network:
        // this compares topologies, whatever quartets the parameters were estimated from
        Network e = est.deepcopynetwork();
        semidirectnetwork(e);
        for (EdgeId x : e.edge)
            if (!e.N(e.getchild(x)).leaf && !(e.E(x).length >= 0.0)) e.E(x).length = 0.0;
        double s_refit = fitnumericalparameters_all(e, q, 0.0, fo, nthreads);
        s_refit = std::max(s_refit, s_est);

        std::printf("{\"hwcd\": %d, \"score_estimated\": %s, \"score_truth\": %s, \"score_refit\": %s, "
                    "\"nhybrids\": %d}\n",
                    hwcd, julia_float(s_est).c_str(), julia_float(s_truth).c_str(), julia_float(s_refit).c_str(),
                    est.numhybrids);
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
