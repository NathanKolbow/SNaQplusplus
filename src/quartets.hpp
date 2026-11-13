// Quartet ranking/sampling and observed concordance factors from gene trees.
#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include "network.hpp"

namespace csnaq {

// Quartets are ranked 0..C(n,4)-1 in colexicographic order of their sorted
// 0-based taxon indices (SNaQ's `unrank4taxa`/`whichQuartet`, shifted to 0).
std::array<int, 4> unrank4taxa(int ntaxa, long long rank);
long long rank4taxa(const std::array<int, 4>& sortedidx);

// SNaQ's `quartetsamplesize(n, p)`.
long long quartetsamplesize(long long n, double p);
// SNaQ's `sampleqindices(n, p, rng)`: sorted 0-based ranks.
std::vector<long long> sampleqindices(long long n, double p, Rng& rng);
// SNaQ's `sampleqindices(n, p, informative, rng)`.
std::vector<long long> sampleqindices(long long n, double p, const std::vector<char>& informative, Rng& rng);

// SNaQ's `isquartetinformative(row, atol)`.
bool isquartetinformative(const double* row, double atol);

// Observed CFs computed from gene trees via a table of pairwise LCA depths
// (SNaQ's `LazyQuartetCF`). Rows are CFs (12|34, 13|24, 14|23) of the quartet
// of the given rank among the sorted taxa.
class ObservedCFs {
public:
    // `lazy_semantics`: if true, a quartet unresolved in a gene tree (polytomy)
    // is counted for the first tied topology, as SNaQ's lazy DataCF does.
    // If false, unresolved quartets are skipped, as the non-lazy DataCF does.
    ObservedCFs(const std::vector<Network>& trees, std::vector<std::string> taxa, bool lazy_semantics);

    const std::vector<std::string>& taxa() const { return taxa_; }
    int ntaxa() const { return static_cast<int>(taxa_.size()); }
    int ntrees() const { return ntrees_; }
    long long nquartets() const { return nchoose4(ntaxa()); }

    // CFs of the quartet with this rank.
    std::array<double, 3> compute(long long rank) const;
    // Fill `out` (3 doubles per rank) for the given ranks; multithreaded.
    void compute_rows(const std::vector<long long>& ranks, double* out, int nthreads) const;
    // Proportion of the gene trees that are counted in the CFs of the quartet
    // with this rank (they have its four taxa and, unless lazy, resolve it).
    double genefraction(long long rank) const;

private:
    std::vector<std::string> taxa_;
    int ntrees_ = 0;
    bool lazy_ = false;
    // lcadepth_[pair * ntrees + t]: depth of the LCA of the pair in tree t,
    // -1 if a taxon of the pair is absent from that tree.
    std::vector<int16_t> lcadepth_;
    long long pairindex(int i, int j) const;  // i < j, 0-based
    void counts(long long rank, long long n[3]) const;
};

// All taxon names appearing in a set of trees, sorted (Julia `sort`).
std::vector<std::string> sorted_union_taxa(const std::vector<Network>& trees);

}  // namespace csnaq
