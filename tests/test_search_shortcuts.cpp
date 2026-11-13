// Tests of the search shortcuts that are not in SNaQ.jl: the score upper
// bound (boundPrune) and the topology signature (dedup).
#include <cmath>

#include "classes.hpp"
#include "signature.hpp"
#include "test_support.hpp"

using namespace csnaq;
using namespace csnaq::test;

namespace {

CFMatrix noisy_cfs(const Network& truth, Rng& rng, double noise) {
    Network c = truth.clone();
    CFMatrix q = computeexpectedCFmatrix(c);
    for (size_t j = 0; j < q.nrow(); ++j) {
        double s = 0.0;
        for (int k = 0; k < 3; ++k) {
            q(j, k) = std::max(1e-3, q(j, k) + noise * (rng.rand() - 0.5));
            s += q(j, k);
        }
        for (int k = 0; k < 3; ++k) q(j, k) /= s;
    }
    return q;
}

Network semidirected(const std::string& newick) {
    Network n = nw(newick).deepcopynetwork();
    semidirectnetwork(n);
    return n;
}

uint64_t sig(const std::string& newick, bool* discrete = nullptr) {
    Network n = semidirected(newick);
    bool d = false;
    uint64_t s = topologysignature(n, &d);
    if (discrete) *discrete = d;
    return s;
}

}  // namespace

TEST(BoundPrune, BoundIsNeverBelowAFittedScore) {
    Rng rng(5);
    Network truth = read_trees("gennet_10_2_8.net").at(0);
    CFMatrix q = noisy_cfs(truth, rng, 0.2);
    const double qlogq = obsqlogq(q);
    Network N = random_tree(10, rng).deepcopynetwork();
    {  // same taxa as the data
        std::vector<std::string> names = truth.tiplabels();
        std::sort(names.begin(), names.end());
        std::vector<NodeId> leaves = N.leaf;
        for (size_t k = 0; k < leaves.size(); ++k) N.N(leaves[k]).name = names[k];
    }
    semidirectnetwork(N);
    N.compact_pids();
    std::vector<AttemptedMove> none;
    std::vector<double> deltas;
    CFMatrix emptyQ;
    int checked = 0, tight = 0;
    for (int it = 0; it < 150; ++it) {
        Network P = N.deepcopynetwork();
        try {
            auto mp = generatemoveproposal(P, nullptr, none, 2, 0.0, emptyQ, deltas, rng, 0.0);
            if (mp.first == Move::none) continue;
            applymove(P, mp.first, mp.second, rng);
            removedegree2nodes(P);
            while (shrink3cycles(P)) {
            }
            while (shrink2cycles(P)) {
            }
            semidirectnetwork(P);
        } catch (const std::exception&) {
            continue;
        }
        for (EdgeId e : P.edge)
            if (!P.N(P.getchild(e)).leaf) P.E(e).length = std::max(P.E(e).length, 0.0);
        ClassOptions co;
        co.qlogq = qlogq;
        co.defer = true;
        QuartetClasses cls = buildclasses(P, QuartetList::allquartets(P.numtaxa), q, co);
        const double early = classupperbound(cls);  // before the reticulate equations are built
        finishclasses(P, cls);
        const double bound = classupperbound(cls);
        FitOptions fo;
        fo.maxeval = 2000;
        const double fitted = fitclasses(P, cls, 0.0, fo);
        ASSERT_GE(early, bound - 1e-9);
        ASSERT_GE(bound, fitted - 1e-9) << "iteration " << it;
        ++checked;
        if (bound - fitted < 0.5 * std::abs(fitted)) ++tight;
        if (rng.rand() < 0.5) N = std::move(P);
    }
    EXPECT_GT(checked, 100);
    EXPECT_GT(tight, 0);
}

TEST(BoundPrune, DoesNotChangeTheSearch) {
    Rng rng(9);
    Network truth = read_trees("gennet_10_2_8.net").at(0);
    CFMatrix q = noisy_cfs(truth, rng, 0.1);
    std::vector<std::string> taxa = truth.tiplabels();
    std::sort(taxa.begin(), taxa.end());
    Network start = truth.clone();
    semidirectnetwork(start);
    removeexcessreticulations(start, 0);
    auto run = [&](bool prune) {
        SnaqOptions o = quiet_options();
        o.threads = 1;
        o.hmax = 2;
        o.runs = 3;
        o.Nfail = 15;
        o.dedup = false;
        o.boundPrune = prune;
        CFSource src = CFSource::full(q, taxa);
        std::vector<Network> nets = multisearch({start}, true, src, o.hmax, o, o.propQuartetsFinal);
        std::vector<double> scores;
        for (const Network& n : nets) scores.push_back(n.fscore);
        return scores;
    };
    std::vector<double> with = run(true), without = run(false);
    ASSERT_EQ(with.size(), without.size());
    for (size_t k = 0; k < with.size(); ++k) EXPECT_EQ(with[k], without[k]);
}

TEST(Dedup, SignatureIgnoresRootAndRepresentation) {
    bool discrete = false;
    // one tree, written from three different roots and with children in different orders
    const uint64_t a = sig("((a,b),(c,d),((e,f),(g,h)));", &discrete);
    EXPECT_TRUE(discrete);
    EXPECT_EQ(a, sig("(((h,g),(f,e)),(d,c),(b,a));"));
    EXPECT_EQ(a, sig("(a,b,((c,d),((e,f),(g,h))));"));
    EXPECT_EQ(a, sig("((g,h),(e,f),((a,b),(c,d)));"));
    // different trees
    EXPECT_NE(a, sig("((a,c),(b,d),((e,f),(g,h)));"));
    EXPECT_NE(a, sig("((a,b),(c,e),((d,f),(g,h)));"));
    // one network from two roots; which hybrid edge is major and the γ values do not matter
    const uint64_t n1 = sig("((a,(b)#H1:::0.8),(c,#H1:::0.2),((e,f),(g,h)));", &discrete);
    EXPECT_TRUE(discrete);
    EXPECT_EQ(n1, sig("((a,(b)#H1:::0.3),(c,#H1:::0.7),((e,f),(g,h)));"));
    EXPECT_EQ(n1, sig("((e,f),(g,h),((a,(b)#H1:::0.8),(c,#H1:::0.2)));"));
    EXPECT_NE(n1, a);
    // the direction of the reticulation matters
    EXPECT_NE(n1, sig("(((a)#H1:::0.8,b),(c,#H1:::0.2),((e,f),(g,h)));"));
    // a reticulation between other edges
    EXPECT_NE(n1, sig("((a,(b)#H1:::0.8),c,((e,f),((g,#H1:::0.2),h)));"));
}

TEST(Dedup, SignatureSeparatesTheNeighborsOfATree) {
    // every valid rNNI(1) of a tree gives a different topology, or the tree itself
    Rng rng(3);
    Network N = random_tree(9, rng).deepcopynetwork();
    semidirectnetwork(N);
    bool discrete = false;
    const uint64_t s0 = topologysignature(N, &discrete);
    ASSERT_TRUE(discrete);
    std::vector<uint64_t> sigs;
    for (const auto& c : allvalidrNNInodes(N, 1)) {
        Network P = N.clone();
        try {
            performrNNI(P, c[0], c[1], c[2], c[3], 1, rng);
            semidirectnetwork(P);
        } catch (const std::exception&) {
            continue;
        }
        sigs.push_back(topologysignature(P, &discrete));
        EXPECT_TRUE(discrete);
        // equal signatures must mean equal topologies
        EXPECT_EQ(sigs.back() == s0, hardwiredclusterdistance(P, N, false) == 0);
    }
    std::sort(sigs.begin(), sigs.end());
    const size_t distinct = static_cast<size_t>(std::unique(sigs.begin(), sigs.end()) - sigs.begin());
    EXPECT_EQ(distinct, 2 * (9 - 3));  // an unrooted binary tree on n taxa has 2(n-3) NNI neighbors
}

TEST(Dedup, SearchStillFindsTheNetwork) {
    Rng rng(9);
    Network truth = read_trees("gennet_10_2_8.net").at(0);
    CFMatrix q = noisy_cfs(truth, rng, 0.0);
    std::vector<std::string> taxa = truth.tiplabels();
    std::sort(taxa.begin(), taxa.end());
    Network start = truth.clone();
    semidirectnetwork(start);
    removeexcessreticulations(start, 0);
    for (bool dedup : {false, true}) {
        SnaqOptions o = quiet_options();
        o.threads = 2;
        o.hmax = 2;
        o.runs = 6;
        o.Nfail = 25;
        o.dedup = dedup;
        Network best = snaq_cfs({start}, true, q, taxa, o);
        EXPECT_GT(best.fscore, -0.05) << "dedup = " << dedup;
    }
}
