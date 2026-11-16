// Validation tool: for network(s) read from a Newick file, compute expected
// quartet CFs, SNaQ scores and gradients (at ρ = 0 and ρ = 0.5) against
// synthetic observed CFs, and write them in a binary format compared against
// SNaQ.jl by tools/julia/check_networks.jl.
//
// usage: snaqpp_check <networks file> <line (1-based)> <output file> [max quartets]
//
// Binary layout (little endian):
//   int64 nq, int64 nparam,
//   float64 score_rho0, float64 score_rho0_batch, float64 score_rho05,
//   nparam x float64 grad_rho0, nparam x float64 grad_rho0_batch, nparam x float64 grad_rho05,
//   nq x 3 float64 eCF (ρ = 0), nq x 3 float64 eCF (ρ = 0.5)
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "equations.hpp"
#include "likelihood.hpp"
#include "newick.hpp"
#include "quartets.hpp"
#include "snaq_network.hpp"

using namespace snaqpp;

// Synthetic observed CFs: must match tools/julia/check_networks.jl
static void synthetic_obs(long long r, double* row) {
    double a = 1 + ((r * 7) % 11), b = 1 + ((r * 13 + 5) % 7), c = 1 + ((r * 3 + 1) % 5);
    if (r % 4 == 0) c = 0;
    if (r % 9 == 0) a = 0;
    double s = a + b + c;
    row[0] = a / s;
    row[1] = b / s;
    row[2] = c / s;
}

// Deterministic subset of ranks used when C(n,4) exceeds the cap
std::vector<long long> quartet_subset(long long nq, long long cap) {
    std::vector<long long> r;
    if (cap <= 0 || nq <= cap) {
        r.resize(static_cast<size_t>(nq));
        for (long long i = 0; i < nq; ++i) r[static_cast<size_t>(i)] = i;
    } else {
        r.resize(static_cast<size_t>(cap));
        for (long long i = 0; i < cap; ++i) r[static_cast<size_t>(i)] = (i * nq) / cap;
    }
    return r;
}

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cerr << "usage: snaqpp_check <networks file> <line> <output file> [max quartets] [threads]\n";
        return 2;
    }
    try {
        std::string path = argv[1];
        long line = std::stol(argv[2]);
        std::string outpath = argv[3];
        long long cap = argc > 4 ? std::stoll(argv[4]) : 0;
        int nthreads = argc > 5 ? std::stoi(argv[5]) : 1;
        std::ifstream in(path);
        std::string s;
        for (long k = 0; k < line && std::getline(in, s); ++k) {
        }
        Network net0 = readnewick_string(s);
        Network net = net0.deepcopynetwork();
        semidirectnetwork(net);
        long long nqall = nchoose4(net.numtaxa);
        std::vector<long long> ranks = quartet_subset(nqall, cap);
        EquationSet es = findquartetequations(net, ranks, nthreads);
        const size_t nq = ranks.size();
        CFMatrix obs;
        obs.v.resize(3 * nq);
        for (size_t j = 0; j < nq; ++j) synthetic_obs(ranks[j], &obs.v[3 * j]);
        const std::vector<double>& params = es.info.params;
        const size_t np = params.size();

        std::vector<double> g0, g05;
        double s0 = computelossandgradient(es.eqns, params, g0, obs, rhotoalpha(0.0), nthreads);
        double s05 = computelossandgradient(es.eqns, params, g05, obs, rhotoalpha(0.5), nthreads);
        TreeQuartetBatch batch = treequartetbatch(es.eqns, obs);
        LossWorkspace ws;
        std::vector<double> gb(np);
        double sb = computelossandgradient(batch, es.eqns, params.data(), np, gb.data(), obs, rhotoalpha(0.0), ws,
                                           nthreads);

        std::ofstream out(outpath, std::ios::binary);
        auto w64 = [&](int64_t v) { out.write(reinterpret_cast<const char*>(&v), 8); };
        auto wd = [&](double v) { out.write(reinterpret_cast<const char*>(&v), 8); };
        w64(static_cast<int64_t>(nq));
        w64(static_cast<int64_t>(np));
        wd(s0);
        wd(sb);
        wd(s05);
        for (double v : g0) wd(v);
        for (double v : gb) wd(v);
        for (double v : g05) wd(v);
        for (double rho : {0.0, 0.5}) {
            for (size_t j = 0; j < nq; ++j) {
                auto e = computeexpectedCF(*es.eqns[j], params, rho);
                wd(e[0]);
                wd(e[1]);
                wd(1.0 - e[0] - e[1]);
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
