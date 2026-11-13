// Hill-climbing search over network topologies (topology_search.jl) and
// multiple independent runs (multisearch.jl).
#pragma once

#include <cstdint>
#include <mutex>
#include <string>
#include <unordered_map>
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
    // `weights`: one weight per quartet, or empty (all quartets weigh 1).
    static CFSource full(CFMatrix m, std::vector<std::string> taxa, std::vector<double> weights = {});
    static CFSource lazy(const ObservedCFs* obs);

    bool is_lazy() const { return lazy_ != nullptr; }
    long long nrows() const;
    const std::vector<std::string>& taxa() const;
    // Rows for the given 0-based ranks.
    CFMatrix rows(const std::vector<long long>& ranks, int nthreads) const;
    // The full matrix (non-lazy only).
    const CFMatrix& matrix() const;
    // Weights of the quartets with the given ranks: the proportion of gene
    // trees counted in their CFs. Empty if the source has no weights.
    std::vector<double> weights(const std::vector<long long>& ranks) const;

private:
    CFMatrix full_;
    std::vector<double> fullweights_;
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
    int opt_maxeval = -1;  // -1: max(30, number of taxa) with lbfgs, 200 with newton
    std::string logprefix;
    // CSNaQ options (not in SNaQ.jl)
    // Reject a proposed network without fitting it when an upper bound of its
    // score shows that it cannot be accepted. The proposal is still counted
    // as a failed proposal, so this does not change what a search does.
    bool boundPrune = true;
    // Do not fit a proposed network whose topology this run has already
    // fitted (such proposals are drawn again, like moves that were already
    // attempted). When every run uses the same quartets (propQuartets = 1),
    // topologies fitted by other runs are not fitted again either.
    bool dedup = true;
    // Whether a proposal skipped by `dedup` counts as a failed proposal
    // (towards Nfail). If false, another proposal is drawn instead, so a run
    // tries Nfail distinct topologies before it stops.
    bool dedupCounts = false;
    // Half of the proposed hybrid additions are drawn between edges where the
    // data suggest a reticulation (see hybridcandidates); the other half
    // uniformly, as SNaQ.jl does.
    bool guideHybrids = true;
    // A run whose network has fewer than hmax reticulations does not stop
    // before this many hybrid additions have been evaluated since its last
    // accepted move. -1: Nfail. 0: no such requirement, as SNaQ.jl.
    int minAddHybrid = -1;
    // Weigh each quartet by the proportion of gene trees that its observed
    // CFs come from (1 for all quartets when every gene tree has every taxon),
    // instead of weighing all quartets equally as SNaQ.jl does.
    bool geneWeights = false;
    // How branch lengths and γ are optimized: "newton" (projected Newton on
    // the quartet classes, run to convergence) or "lbfgs" (NLopt's L-BFGS
    // limited to opt_maxeval evaluations, as SNaQ.jl).
    std::string optimizer = "newton";
    // Shallow evaluation while a proposed network is fitted (0: exact; see
    // FitOptions::shallow). `shallowMode`: "path", "average" or "freeze". The
    // score of every network is still its exact score, and the final
    // optimization of a run's network is exact.
    double shallow = 0.0;
    std::string shallowMode = "path";
    // CLI-only
    bool lazy = false;
    int threads = 1;
};

// Networks are compared by how close their score is to 0, the score of a
// model whose expected CFs equal the observed CFs: larger is better. The score
// is at most 0 (up to rounding: expected CFs are floored at 1e-9), so this is
// the order of the scores themselves unless a score is positive.
inline double scorerank(double score) { return score > 0.0 ? -score : score; }

// Scores of the topologies fitted by the runs of a multisearch, shared when
// every run uses the same quartets.
struct SearchShared {
    std::mutex m;
    std::unordered_map<uint64_t, double> score;  // topology signature -> best converged score
};

// One search (SNaQ's `search`). `logfile` and `filename` as in SNaQ.
Network search(Network N, const CFSource& q, int hmax, const SnaqOptions& o, long long seed,
               double propQuartetsFinal, const std::string& logfile, const std::string& filename, int nthreads,
               SearchShared* shared = nullptr);

// SNaQ's `multisearch`: returns the networks of all runs, sorted by score
// (best first). `startnets` has length 1 or `o.runs`.
std::vector<Network> multisearch(std::vector<Network> startnets, bool single_input, const CFSource& q, int hmax,
                                 const SnaqOptions& o, double propQuartetsFinal);

// rescoreonsharedquartets!(nets, q, propQuartetsFinal, seed, ρ)
long long rescoreonsharedquartets(std::vector<Network>& nets, const CFSource& q, double propQuartetsFinal,
                                  long long seed, double rho, int nthreads, bool geneweights = false);

// Append `msg` and a newline to "<filename>.log" (SNaQ's logmessage).
void logmessage(const std::string& filename, const std::string& msg);

}  // namespace csnaq
