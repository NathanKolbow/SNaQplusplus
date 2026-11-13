// Hill-climbing search over network topologies (topology_search.jl) and
// multiple independent runs (multisearch.jl).
#pragma once

#include <string>
#include <vector>

#include "likelihood.hpp"
#include "network.hpp"
#include "quartets.hpp"
#include "restrictions.hpp"

namespace csnaq {

// Observed CFs of all C(n,4) quartets, either materialized (non-lazy DataCF)
// or computed on demand from gene trees (lazy DataCF).
class CFSource {
public:
    static CFSource full(CFMatrix m, std::vector<std::string> taxa);
    static CFSource lazy(const ObservedCFs* obs);

    bool is_lazy() const { return lazy_ != nullptr; }
    long long nrows() const;
    const std::vector<std::string>& taxa() const;
    // Rows for the given 0-based ranks.
    CFMatrix rows(const std::vector<long long>& ranks, int nthreads) const;
    // The full matrix (non-lazy only).
    const CFMatrix& matrix() const;

private:
    CFMatrix full_;
    std::vector<std::string> taxa_;
    const ObservedCFs* lazy_ = nullptr;
};

// Every keyword argument of `snaq!` plus the search keyword arguments it
// passes through (`kwargs...`).
struct SnaqOptions {
    int hmax = 1;
    int Nfail = 50;
    double ftolRel = 1e-8;
    double ftolAbs = 1e-8;
    double xtolRel = 1e-8;
    double xtolAbs = 1e-8;
    bool verbose = false;
    int runs = 100;
    std::string outgroup = "none";
    std::string filename = "snaq";
    long long seed = 0;
    bool seed_given = false;
    double probST = 0.9;
    bool updateBL = true;
    double probQR = 0.0;
    double qtolAbs = 1e-4;
    bool qinfTest = false;
    double propQuartets = 1.0;
    double propQuartetsFinal = 1.0;
    Restrictions restrictions;
    double rho = 0.0;
    // search keyword arguments reachable through snaq!'s kwargs...
    long long maxeval = 100000000;
    double liktolAbs = 1e-8;
    double liktolRel = 1e-4;
    int opt_maxeval = -1;  // -1: max(30, number of taxa)
    std::string logprefix;
    // CLI-only
    bool lazy = false;
    int threads = 1;
};

// One search (SNaQ's `search`). `logfile` and `filename` as in SNaQ.
Network search(Network N, const CFSource& q, int hmax, const SnaqOptions& o, long long seed,
               double propQuartetsFinal, const std::string& logfile, const std::string& filename, int nthreads);

// SNaQ's `multisearch`: returns the networks of all runs, sorted by score
// (best first). `startnets` has length 1 or `o.runs`.
std::vector<Network> multisearch(std::vector<Network> startnets, bool single_input, const CFSource& q, int hmax,
                                 const SnaqOptions& o, double propQuartetsFinal);

// rescoreonsharedquartets!(nets, q, propQuartetsFinal, seed, ρ)
long long rescoreonsharedquartets(std::vector<Network>& nets, const CFSource& q, double propQuartetsFinal,
                                  long long seed, double rho, int nthreads);

// Append `msg` and a newline to "<filename>.log" (SNaQ's logmessage).
void logmessage(const std::string& filename, const std::string& msg);

}  // namespace csnaq
