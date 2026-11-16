#include "snaq.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <random>
#include <set>

#include "classes.hpp"
#include "equations.hpp"
#include "network_manip.hpp"
#include "optimize.hpp"
#include "quartets.hpp"
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
    // As lazysnaq!: no parameter re-optimization after a run (propQuartetsFinal
    // = 0 for the runs); the networks of all runs are compared on one shared
    // sample of the quartets, with the parameters fitted during their search.
    // (SNaQ.jl sorts the output files by the per-run scores, which are not
    // comparable, and only returns the network selected on the shared quartets;
    // here the output files are sorted by the scores on the shared quartets,
    // so their first network is the one that is returned.)
    std::vector<Network> all_nets = multisearch(startnets, single_input, q, o.hmax, o, 0.0, o.propQuartetsFinal);
    logmessage(o.filename, "snaq! (lazy DataCF): selected the best network after fair comparison of the " +
                               std::to_string(o.runs) + " runs on shared quartets (SNaQscore = " +
                               julia_float(julia_round(all_nets[0].fscore, 5)) + ").");
    return finish(std::move(all_nets[0]));
}

}  // namespace

static void default_seed(SnaqOptions& o) {
    if (!o.seed_given) {
        // snaq!: seed::Int = rand(Int)
        std::random_device rd;
        uint64_t s = (static_cast<uint64_t>(rd()) << 32) ^ rd();
        o.seed = static_cast<long long>(s);
        o.seed_given = true;
    }
}

Network snaq_cfs(const std::vector<Network>& startnets, bool single_input, CFMatrix cfs,
                 std::vector<std::string> taxa, SnaqOptions o) {
    if (startnets.empty()) throw SnaqError("no starting topology given");
    default_seed(o);
    o.threads = std::max(1, o.threads);
    std::sort(taxa.begin(), taxa.end());
    if (static_cast<long long>(cfs.nrow()) != nchoose4(static_cast<long long>(taxa.size())))
        throw SnaqError("the CF matrix must have one row per quartet of the taxa");
    CFSource q = CFSource::full(std::move(cfs), taxa);
    std::vector<Network> sorted = multisearch(startnets, single_input, q, o.hmax, o, o.propQuartetsFinal);
    return finish(std::move(sorted[0]));
}

// Whether lazysnaq accepts these inputs: there are gene trees, and every
// starting topology has exactly their taxa.
static bool lazyapplies(const std::vector<Network>& startnets, const std::vector<Network>& genetrees) {
    if (genetrees.empty()) return false;
    const std::vector<std::string> taxa = sorted_union_taxa(genetrees);
    for (const Network& n : startnets) {
        std::vector<std::string> t = n.tiplabels();
        std::sort(t.begin(), t.end());
        if (t != taxa) return false;
    }
    return true;
}

Network snaq(const std::vector<Network>& startnets, bool single_input, const std::vector<Network>& genetrees,
             SnaqOptions o) {
    if (startnets.empty()) throw SnaqError("no starting topology given");
    default_seed(o);
    o.threads = std::max(1, o.threads);
    if (o.lazy) {
        if (o.propQuartets == 1.0)
            throw SnaqError("snaq! with a lazy DataCF requires propQuartets < 1.0, otherwise every quartet is "
                            "computed anyways. Specify propQuartets and propQuartetsFinal, or use a DataCF with "
                            "lazy=false.");
        return lazysnaq(startnets, single_input, genetrees, o);
    }
    // Not asked for, but only a sample of the quartets is ever needed: observed
    // CFs are computed on demand, where the lazy DataCF can be used.
    if (o.autoLazy && 0 < o.propQuartets && o.propQuartets < 1 && 0 < o.propQuartetsFinal &&
        o.propQuartetsFinal < 1 && !o.qinfTest && lazyapplies(startnets, genetrees)) {
        o.lazy = true;
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
    std::vector<double> weights;
    if (o.geneWeights) {
        weights.resize(static_cast<size_t>(nq));
        for (long long r = 0; r < nq; ++r) weights[static_cast<size_t>(r)] = obs.genefraction(r);
    }
    CFSource q = CFSource::full(std::move(m), obs.taxa(), std::move(weights));

    std::vector<Network> sorted = multisearch(startnets, single_input, q, o.hmax, o, o.propQuartetsFinal);
    return finish(std::move(sorted[0]));
}

Network optimize_network(Network net, const std::vector<Network>& genetrees, SnaqOptions o) {
    default_seed(o);
    const int nthreads = std::max(1, o.threads);
    if (!(0 < o.propQuartets && o.propQuartets <= 1))
        throw SnaqError("propQuartets must be in range (0, 1] (propQuartets = " + julia_float(o.propQuartets) + ")");
    if (!(0 <= o.rho && o.rho <= 1)) throw SnaqError("ρ must be in range [0, 1] (ρ = " + julia_float(o.rho) + ")");
    FitOptions fit;
    if (o.optimizer == "newton")
        fit.method = FitMethod::newton;
    else if (o.optimizer == "lbfgs")
        fit.method = FitMethod::lbfgs;
    else
        throw SnaqError("optimizer must be \"newton\" or \"lbfgs\" (optimizer = " + o.optimizer + ")");
    fit.maxeval = o.opt_maxeval;
    fit.ftolRel = o.ftolRel;
    fit.ftolAbs = o.ftolAbs;
    fit.xtolRel = o.xtolRel;
    fit.xtolAbs = o.xtolAbs;
    if (genetrees.empty()) throw SnaqError("no gene trees given");

    // taxa: those of the gene trees; the network must have all of them
    const std::vector<std::string> taxa = sorted_union_taxa(genetrees);
    net = net.deepcopynetwork();
    {
        auto tl = net.tiplabels();
        std::set<std::string> tset(tl.begin(), tl.end()), dset(taxa.begin(), taxa.end());
        std::string missing, extra;
        for (const auto& t : taxa)
            if (!tset.count(t)) missing += (missing.empty() ? "" : ", ") + ("\"" + t + "\"");
        if (!missing.empty()) throw SnaqError("Gene trees contain taxa not present in the network: [" + missing + "].");
        for (const auto& t : tl)
            if (!dset.count(t)) extra += (extra.empty() ? "" : ", ") + ("\"" + t + "\"");
        if (!extra.empty()) {
            warn("The following taxa are in the network but in no gene tree, they will be deleted: [" + extra + "]");
            for (const auto& t : tl)
                if (!dset.count(t)) deleteleaf_name(net, t, DeleteOpts{});
        }
    }
    for (NodeId n : net.node)
        if (net.N(n).edge.size() > 3) throw SnaqError("the network has a polytomy: SNaQ only fits binary networks");
    if (net.numtaxa < 4) throw SnaqError("the network must have at least 4 taxa");

    // as for the starting topologies of snaq!: missing lengths are 0, invalid γ are 0.5
    semidirectnetwork(net);
    for (EdgeId e : net.edge)
        if (net.E(e).length == -1.0) net.E(e).length = 0.0;
    for (NodeId h : net.hybrid) {
        double gM = net.E(net.getparentedge(h)).gamma, gm = net.E(net.getparentedgeminor(h)).gamma;
        bool sum1 = std::abs(gM + gm - 1.0) <= 1.4901161193847656e-8 * std::max(std::abs(gM + gm), 1.0);
        if (1 >= gM && gM >= 0 && 1 >= gm && gm >= 0 && sum1) continue;
        net.E(net.getparentedge(h)).gamma = 0.5;
        net.E(net.getparentedgeminor(h)).gamma = 0.5;
    }
    for (NodeId n : net.node)
        if (!net.N(n).leaf && !net.N(n).hybrid) net.N(n).name.clear();

    // observed CFs of the used quartets only
    ObservedCFs obs(genetrees, taxa, /*lazy_semantics=*/false);
    const long long nq = obs.nquartets();
    std::vector<long long> ranks;
    if (o.propQuartets == 1.0) {
        ranks.resize(static_cast<size_t>(nq));
        std::iota(ranks.begin(), ranks.end(), 0LL);
    } else {
        Rng rng(static_cast<uint64_t>(o.seed));
        ranks = sampleqindices(nq, o.propQuartets, rng);
    }
    std::cout << "calculating obsCF from " << genetrees.size() << " gene trees and for " << ranks.size() << " of "
              << nq << " quartets" << std::endl;
    CFMatrix m;
    m.v.resize(3 * ranks.size());
    obs.compute_rows(ranks, m.v.data(), nthreads);
    std::vector<double> weights;
    if (o.geneWeights) {
        weights.resize(ranks.size());
        for (size_t k = 0; k < ranks.size(); ++k) weights[k] = obs.genefraction(ranks[k]);
    }

    ClassOptions co;
    co.weights = weights.empty() ? nullptr : weights.data();
    co.qlogq = obsqlogq(m, co.weights);
    co.nthreads = nthreads;
    fixnegativeedges(net);
    checkgammas(net);
    net.compact_pids();
    QuartetClasses cls = buildclasses(net, QuartetList::fromranks(net.numtaxa, ranks), m, co);
    const double score = fitclasses(net, cls, o.rho, fit, nthreads);

    for (NodeId n : net.node)
        if (!net.N(n).leaf && !net.N(n).hybrid) net.N(n).name.clear();
    for (size_t iH = 0; iH < net.hybrid.size(); ++iH) net.N(net.hybrid[iH]).name = "H" + std::to_string(iH + 1);
    if (o.outgroup != "none") {
        bool has = false;
        for (NodeId l : net.leaf) has |= net.N(l).name == o.outgroup;
        if (!has) throw ArgumentError("The network does not contain the supplied outgroup (" + o.outgroup + ").");
        try {
            rootatnode_name(net, o.outgroup);
        } catch (const RootMismatch&) {
            throw ArgumentError("The network contains the outgroup but cannot be rooted at the outgroup.");
        }
    }
    net = finish(std::move(net));
    net.fscore = score;
    return net;
}

}  // namespace csnaq
