#include "snaq.hpp"

#include <algorithm>
#include <iostream>
#include <numeric>
#include <random>
#include <set>

#include "snaq_network.hpp"

namespace csnaq {

namespace {

Network finish(Network bestnet) {
    semidirectnetwork(bestnet);
    for (NodeId l : bestnet.leaf) bestnet.E(bestnet.getparentedge(l)).length = 0.0;
    return bestnet;
}

// lazysnaq!(currT0, d, hmax, propQuartets, propQuartetsFinal; ...)
Network lazysnaq(const std::vector<Network>& startnets, bool single_input, const std::vector<Network>& genetrees,
                 const SnaqOptions& o) {
    if (o.qinfTest)
        throw SnaqError("snaq! with a lazy DataCF does not support qinfTest=true: it needs every quartet's observed "
                        "CF up front. Use a DataCF with lazy=false instead.");
    if (!(0 < o.propQuartets && o.propQuartets <= 1))
        throw SnaqError("propQuartets must be in range (0, 1] (propQuartets = " + julia_float(o.propQuartets) + ")");
    if (!(0 < o.propQuartetsFinal && o.propQuartetsFinal <= 1))
        throw SnaqError("propQuartetsFinal must be in range (0, 1] (propQuartetsFinal = " +
                        julia_float(o.propQuartetsFinal) + ")");
    if (genetrees.empty())
        throw SnaqError("snaq! with a lazy DataCF needs the gene trees its CFs are computed from, but it has none.");
    if (o.propQuartetsFinal == 1.0)
        warn("propQuartetsFinal is 1.0 (the default value). This may be very time consuming for large datasets.");

    std::vector<std::string> taxa = startnets[0].tiplabels();
    std::sort(taxa.begin(), taxa.end());
    if (!single_input)
        for (const Network& n : startnets) {
            auto t = n.tiplabels();
            std::sort(t.begin(), t.end());
            if (t != taxa)
                throw SnaqError("All starting topologies in currT0 must share the same set of taxa, so that one "
                                "LazyQuartetCF's rank <-> taxa mapping is valid for every run.");
        }
    // DataCF(genetrees; lazy=true): taxa are those of the gene trees
    ObservedCFs lazyq(genetrees, sorted_union_taxa(genetrees), /*lazy_semantics=*/true);
    std::set<std::string> tset(taxa.begin(), taxa.end()), lset(lazyq.taxa().begin(), lazyq.taxa().end());
    std::string extra, nevershown;
    for (const auto& t : lazyq.taxa())
        if (!tset.count(t)) extra += (extra.empty() ? "" : ", ") + ("\"" + t + "\"");
    if (!extra.empty()) throw SnaqError("Gene trees contain taxa not present in currT0: [" + extra + "].");
    for (const auto& t : taxa)
        if (!lset.count(t)) nevershown += (nevershown.empty() ? "" : ", ") + ("\"" + t + "\"");
    if (!nevershown.empty())
        throw SnaqError("These taxa in currT0 appear in no gene tree, so every quartet involving them would have an "
                        "observed CF of NaN: [" + nevershown + "].");

    CFSource q = CFSource::lazy(&lazyq);
    std::vector<Network> all_nets = multisearch(startnets, single_input, q, o.hmax, o, 0.0);
    long long nusedfinal = rescoreonsharedquartets(all_nets, q, o.propQuartetsFinal, o.seed, o.rho, o.threads);
    size_t bestidx = 0;
    for (size_t k = 1; k < all_nets.size(); ++k)
        if (all_nets[k].fscore > all_nets[bestidx].fscore) bestidx = k;
    // QUIRK: `bestidx` is the position among the runs' networks sorted by their
    // (not comparable) per-run scores, not the run number.
    logmessage(o.filename, "snaq! (lazy DataCF): selected run " + std::to_string(bestidx + 1) + " of " +
                               std::to_string(o.runs) +
                               " as the best network after fair comparison (SNaQscore = " +
                               julia_float(julia_round(all_nets[bestidx].fscore, 5)) + " on " +
                               std::to_string(nusedfinal) + " shared quartets).");
    return finish(std::move(all_nets[bestidx]));
}

}  // namespace

Network snaq(const std::vector<Network>& startnets, bool single_input, const std::vector<Network>& genetrees,
             SnaqOptions o) {
    if (startnets.empty()) throw SnaqError("no starting topology given");
    if (!o.seed_given) {
        // snaq!: seed::Int = rand(Int)
        std::random_device rd;
        uint64_t s = (static_cast<uint64_t>(rd()) << 32) ^ rd();
        o.seed = static_cast<long long>(s);
        o.seed_given = true;
    }
    o.threads = std::max(1, o.threads);
    if (o.lazy) {
        if (o.propQuartets == 1.0)
            throw SnaqError("snaq! with a lazy DataCF requires propQuartets < 1.0, otherwise every quartet is "
                            "computed anyways. Specify propQuartets and propQuartetsFinal, or use a DataCF with "
                            "lazy=false.");
        return lazysnaq(startnets, single_input, genetrees, o);
    }

    // DataCF(genetrees): observed CFs of every quartet, computed up front.
    std::vector<std::string> taxa = sorted_union_taxa(genetrees);
    ObservedCFs obs(genetrees, taxa, /*lazy_semantics=*/false);
    const long long nq = obs.nquartets();
    std::cout << "calculating obsCF from " << genetrees.size() << " gene trees and for " << nq << " quartets"
              << std::endl;
    std::vector<long long> ranks(static_cast<size_t>(nq));
    std::iota(ranks.begin(), ranks.end(), 0LL);
    CFMatrix m;
    m.v.resize(3 * static_cast<size_t>(nq));
    obs.compute_rows(ranks, m.v.data(), o.threads);
    CFSource q = CFSource::full(std::move(m), obs.taxa());

    std::vector<Network> sorted = multisearch(startnets, single_input, q, o.hmax, o, o.propQuartetsFinal);
    return finish(std::move(sorted[0]));
}

}  // namespace csnaq
