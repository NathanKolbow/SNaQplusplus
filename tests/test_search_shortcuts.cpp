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
        o.runs = 10;
        o.Nfail = 25;
        o.dedup = dedup;
        Network best = snaq_cfs({start}, true, q, taxa, o);
        EXPECT_GT(best.fscore, -0.1) << "dedup = " << dedup;
    }
}

TEST(NewtonFit, ReachesTheOptimumOfLbfgs) {
    Rng rng(31);
    std::vector<Network> nets;
    nets.push_back(nw("(15,(1,((14,(#H1,(((12,13),(11,#H3)),(7,((10)#H3,(8,9)))))),((((2,3))#H2,(6,(5,(#H2,4)))))#H1)));"));
    nets.push_back(read_trees("gennet_10_2_8.net").at(0));
    nets.push_back(read_trees("gennet_10_2_12.net").at(0));
    nets.push_back(random_tree(14, rng));
    for (Network& truth : nets) {
        truth = truth.deepcopynetwork();
        semidirectnetwork(truth);
        for (EdgeId e : truth.edge) truth.E(e).length = 0.2 + rng.rand();
        for (NodeId h : truth.hybrid) {
            truth.E(truth.getparentedgeminor(h)).gamma = 0.3;
            truth.E(truth.getparentedge(h)).gamma = 0.7;
        }
        CFMatrix q = noisy_cfs(truth, rng, 0.1);
        for (int start = 0; start < 3; ++start) {
            Network a = truth.deepcopynetwork(), b = truth.deepcopynetwork();
            for (EdgeId e : a.edge) a.E(e).length = b.E(e).length = start == 0 ? 0.0 : 0.05 + 2.0 * rng.rand();
            ClassOptions co;
            co.qlogq = obsqlogq(q);
            QuartetClasses ca = buildclasses(a, QuartetList::allquartets(a.numtaxa), q, co);
            QuartetClasses cb = buildclasses(b, QuartetList::allquartets(b.numtaxa), q, co);
            FitOptions lb, nt;
            lb.maxeval = 20000;
            nt.method = FitMethod::newton;
            FitInfo ia, ib;
            const double sl = fitclasses(a, ca, 0.0, lb, 1, &ia);
            const double sn = fitclasses(b, cb, 0.0, nt, 1, &ib);
            EXPECT_TRUE(ib.converged);
            // never meaningfully below L-BFGS run to convergence
            EXPECT_GT(sn, sl - 1e-6 * std::max(1.0, std::abs(sl))) << "start " << start;
            if (truth.numhybrids == 0) {  // concave score, exact curvature: few evaluations
                EXPECT_NEAR(sn, sl, 1e-7 * std::max(1.0, std::abs(sl)));
                EXPECT_LT(ib.evals, 30);
            }
        }
    }
}

TEST(NewtonFit, DoesNotStopAtBoundsOrStepLimits) {
    // Optima with branch lengths at their bounds, started far away: Newton
    // steps are cut by the bounds and by the step limits. The fit must still
    // reach the optimum that L-BFGS run to convergence finds, in few evaluations.
    Rng rng(43);
    std::vector<Network> nets;
    nets.push_back(nw("(15,(1,((14,(#H1,(((12,13),(11,#H3)),(7,((10)#H3,(8,9)))))),((((2,3))#H2,(6,(5,(#H2,4)))))#H1)));"));
    nets.push_back(read_trees("gennet_10_2_8.net").at(0));
    nets.push_back(random_tree(16, rng));
    int total = 0, fits = 0;
    for (Network& truth : nets) {
        truth = truth.deepcopynetwork();
        semidirectnetwork(truth);
        for (EdgeId e : truth.edge) {
            const double u = rng.rand();
            truth.E(e).length = u < 0.3 ? 0.0 : (u < 0.6 ? 12.0 : 0.2 + rng.rand());  // many optima at a bound
        }
        for (NodeId h : truth.hybrid) {
            truth.E(truth.getparentedgeminor(h)).gamma = 0.3;
            truth.E(truth.getparentedge(h)).gamma = 0.7;
        }
        CFMatrix q = noisy_cfs(truth, rng, 0.02);
        for (int start = 0; start < 3; ++start) {
            Network a = truth.deepcopynetwork(), b = truth.deepcopynetwork();
            for (EdgeId e : a.edge) a.E(e).length = b.E(e).length = start == 0 ? 1.0 : 0.05 + 3.0 * rng.rand();
            ClassOptions co;
            co.qlogq = obsqlogq(q);
            QuartetClasses ca = buildclasses(a, QuartetList::allquartets(a.numtaxa), q, co);
            QuartetClasses cb = buildclasses(b, QuartetList::allquartets(b.numtaxa), q, co);
            FitOptions lb, nt;
            lb.maxeval = 50000;
            nt.method = FitMethod::newton;
            FitInfo ia, ib;
            const double sl = fitclasses(a, ca, 0.0, lb, 1, &ia);
            const double sn = fitclasses(b, cb, 0.0, nt, 1, &ib);
            EXPECT_TRUE(ib.converged) << "start " << start;
            EXPECT_GT(sn, sl - 1e-6 * std::max(1.0, std::abs(sl))) << "start " << start;
            EXPECT_LT(ib.evals, 100) << "start " << start;
            total += ib.evals;
            ++fits;
        }
    }
    EXPECT_LT(total, 40 * fits);
}

TEST(NewtonFit, ShortensLongBranchesThatTheDataWantShort) {
    // A proposed network inherits its branch lengths from the current one,
    // where some are at the upper bound. Where exp(-length) is close to 0 the
    // derivative with respect to the length is close to 0 too, whatever the
    // score would gain from a shorter branch: the fit must not stop there.
    // Reference: the fit started from the parameters the data were made with.
    Rng rng(47);
    std::vector<Network> nets;
    nets.push_back(nw("(15,(1,((14,(#H1,(((12,13),(11,#H3)),(7,((10)#H3,(8,9)))))),((((2,3))#H2,(6,(5,(#H2,4)))))#H1)));"));
    nets.push_back(read_trees("gennet_10_2_8.net").at(0));
    nets.push_back(random_tree(16, rng));
    int trapped = 0, fits = 0;
    for (Network& truth : nets) {
        truth = truth.deepcopynetwork();
        semidirectnetwork(truth);
        for (EdgeId e : truth.edge) truth.E(e).length = 0.2 + rng.rand();
        for (NodeId h : truth.hybrid) {
            truth.E(truth.getparentedgeminor(h)).gamma = 0.3;
            truth.E(truth.getparentedge(h)).gamma = 0.7;
        }
        CFMatrix q = noisy_cfs(truth, rng, 0.02);
        ClassOptions co;
        co.qlogq = obsqlogq(q);
        FitOptions nt;
        nt.method = FitMethod::newton;
        Network ref = truth.deepcopynetwork();
        QuartetClasses cr = buildclasses(ref, QuartetList::allquartets(ref.numtaxa), q, co);
        const double sref = fitclasses(ref, cr, 0.0, nt, 1);
        for (int start = 0; start < 6; ++start) {
            Network a = truth.deepcopynetwork();
            for (EdgeId e : a.edge)
                if (rng.rand() < 0.3) a.E(e).length = 25.0;
            QuartetClasses ca = buildclasses(a, QuartetList::allquartets(a.numtaxa), q, co);
            FitInfo ia;
            const double sn = fitclasses(a, ca, 0.0, nt, 1, &ia);
            EXPECT_TRUE(ia.converged) << "start " << start;
            if (!(sn > sref - 1e-5 * std::max(1.0, std::abs(sref)))) ++trapped;
            EXPECT_GT(sn, sref - 1e-5 * std::max(1.0, std::abs(sref))) << "start " << start;
            ++fits;
        }
    }
    EXPECT_EQ(trapped, 0) << "of " << fits << " fits";
}

TEST(Sweep, RunStopsWhereNoLocalMoveIsAccepted) {
    // With `sweep`, a run stops only after every local move of its network was
    // proposed without being accepted. Until its first sweep the run is the
    // run without `sweep`, so it cannot end with a worse score.
    Rng rng(21);
    Network truth = read_trees("gennet_10_2_8.net").at(0);
    CFMatrix q = noisy_cfs(truth, rng, 0.1);
    std::vector<std::string> taxa = truth.tiplabels();
    std::sort(taxa.begin(), taxa.end());
    Network start = random_tree(10, rng).deepcopynetwork();
    {  // same taxa as the data
        std::vector<NodeId> leaves = start.leaf;
        for (size_t k = 0; k < leaves.size(); ++k) start.N(leaves[k]).name = taxa[k];
    }
    CFSource src = CFSource::full(q, taxa);
    ClassOptions co;
    co.qlogq = obsqlogq(q);
    int better = 0;
    for (long long seed : {3, 4, 5, 6}) {
        auto run = [&](bool sweep) {
            SnaqOptions o = quiet_options();
            o.hmax = 1;
            o.Nfail = 6;  // few random proposals: runs without the sweep stop early
            o.probST = 0.0;
            o.sweep = sweep;
            return search(start.clone(), src, o.hmax, o, seed, 0.0, "", "", 1);
        };
        Network swept = run(true), plain = run(false);
        EXPECT_GE(swept.fscore, plain.fscore - 1e-12) << "seed " << seed;
        if (swept.fscore > plain.fscore + 1e-6) ++better;
        // no local move of the final network is accepted
        const SnaqOptions o = quiet_options();
        FitOptions fo;
        fo.method = FitMethod::newton;
        fo.ftolRel = o.ftolRel;
        fo.ftolAbs = o.ftolAbs;
        fo.xtolRel = o.xtolRel;
        fo.xtolAbs = o.xtolAbs;
        const double cur = swept.fscore;
        Network N = swept.deepcopynetwork();
        semidirectnetwork(N);
        N.compact_pids();
        size_t nmoves = 0;
        for (const auto& mv : localmoves(N)) {
            Network P = N.deepcopynetwork();
            applymove(P, mv.first, mv.second, rng);
            removedegree2nodes(P);
            while (shrink3cycles(P)) {
            }
            while (shrink2cycles(P)) {
            }
            semidirectnetwork(P);
            QuartetClasses cls = buildclasses(P, QuartetList::allquartets(P.numtaxa), q, co);
            const double score = fitclasses(P, cls, 0.0, fo, 1);
            const bool improves = score - cur > o.liktolAbs && (cur - score) / cur > o.liktolRel;
            EXPECT_FALSE(improves) << "seed " << seed << ", " << move_name(mv.first) << ": " << score << " > " << cur;
            ++nmoves;
        }
        EXPECT_GT(nmoves, 10u);
    }
    EXPECT_GT(better, 0);  // the sweep found something in at least one run
}

TEST(Restrictions, SearchedNetworksMeetThemWhenReadAgain) {
    // The tree-child conditions read attributes of the edges (their direction,
    // whether the root can be on them) that moves leave out of date. A network
    // accepted under a restriction must meet it as it is read from the output:
    // written, read, copied and semidirected.
    Rng rng(31);
    Network truth = nw("(15,(1,((14,(#H1,(((12,13),(11,#H3)),(7,((10)#H3,(8,9)))))),((((2,3))#H2,(6,(5,(#H2,4)))))#H1)));");
    truth = truth.deepcopynetwork();
    semidirectnetwork(truth);
    for (EdgeId e : truth.edge) truth.E(e).length = 0.2 + rng.rand();
    for (NodeId h : truth.hybrid) {
        truth.E(truth.getparentedgeminor(h)).gamma = 0.3;
        truth.E(truth.getparentedge(h)).gamma = 0.7;
    }
    CFMatrix q = noisy_cfs(truth, rng, 0.05);
    std::vector<std::string> taxa = truth.tiplabels();
    std::sort(taxa.begin(), taxa.end());
    Network start = truth.clone();
    semidirectnetwork(start);
    removeexcessreticulations(start, 0);
    CFSource src = CFSource::full(q, taxa);
    size_t checked = 0, withhybrids = 0;
    for (int which = 0; which < 3; ++which) {
        SnaqOptions o = quiet_options();
        o.hmax = 4;
        o.runs = 12;
        o.Nfail = 15;
        o.seed = 5 + which;
        o.restrictions.rooted_tree_child = which == 0;
        o.restrictions.weakly_tree_child = which == 1;
        o.restrictions.strongly_tree_child = which == 2;
        for (Network& net : multisearch({start}, true, src, o.hmax, o, 0.0)) {
            Network again = readnewick_string(writenewick(net)).deepcopynetwork();
            semidirectnetwork(again);
            const auto tc = istreechild(again);
            EXPECT_TRUE(which == 0 ? std::get<0>(tc) : (which == 1 ? std::get<1>(tc) : std::get<2>(tc)))
                << (which == 0 ? "rooted" : (which == 1 ? "weakly" : "strongly")) << " tree-child: " << writenewick(net);
            ++checked;
            if (net.numhybrids >= 2) ++withhybrids;
        }
    }
    EXPECT_EQ(checked, 36u);
    EXPECT_GT(withhybrids, 6u);  // the restrictions had something to restrict
}

TEST(HybridGuide, CandidatesPointAtTheMissingReticulation) {
    // Data from a network with one reticulation, fitted on its major tree: the
    // two edges the reticulation connects must be among the best candidate
    // pairs (the pairs above it connect edges next to them).
    Network truth = nw("((a:1,b:1):1,((c:1,(d:0.5)#H1:0.5::0.7):1,(e:1,(f:1,#H1:0.5::0.3):1):1):1,(g:1,h:1):1);");
    Network t2 = truth.clone();
    CFMatrix q = computeexpectedCFmatrix(t2);
    Network tree = truth.deepcopynetwork();
    semidirectnetwork(tree);
    removeexcessreticulations(tree, 0);
    tree.compact_pids();
    ClassOptions co;
    co.qlogq = obsqlogq(q);
    QuartetClasses cls = buildclasses(tree, QuartetList::allquartets(tree.numtaxa), q, co);
    FitOptions fo;
    fo.method = FitMethod::newton;
    fitclasses(tree, cls, 0.0, fo);
    std::vector<HybridCandidate> cands = hybridcandidates(tree, cls, 0.0);
    ASSERT_FALSE(cands.empty());
    // the reticulation goes between the pendant edges of d and f
    int pd = -1, pf = -1;
    for (NodeId l : tree.leaf) {
        if (tree.N(l).name == "d") pd = tree.E(tree.getparentedge(l)).pid;
        if (tree.N(l).name == "f") pf = tree.E(tree.getparentedge(l)).pid;
    }
    size_t rank = cands.size();
    for (size_t k = 0; k < cands.size(); ++k)
        if ((cands[k].pid1 == pd && cands[k].pid2 == pf) || (cands[k].pid1 == pf && cands[k].pid2 == pd)) rank = k;
    EXPECT_LT(rank, 6u) << "of " << cands.size() << " candidate pairs";
    EXPECT_GT(cands.size(), 20u);
}

TEST(GeneWeights, QuartetsWeighByTheirShareOfGeneTrees) {
    // four gene trees; taxon f is missing from two of them
    std::vector<Network> trees{nw("((a,b),(c,d),(e,f));"), nw("((a,c),(b,d),(e,f));"), nw("((a,b),(c,(d,e)));"),
                               nw("(((a,b),c),(d,e));")};
    std::vector<std::string> taxa{"a", "b", "c", "d", "e", "f"};
    ObservedCFs obs(trees, taxa, false);
    const long long nq = obs.nquartets();
    CFMatrix q;
    q.v.resize(3 * static_cast<size_t>(nq));
    std::vector<long long> ranks;
    std::vector<double> w;
    for (long long r = 0; r < nq; ++r) {
        ranks.push_back(r);
        w.push_back(obs.genefraction(r));
        auto t = unrank4taxa(6, r);
        // quartets with f come from two gene trees (or fewer, if unresolved), the others from up to four
        if (t[3] == 5) {
            EXPECT_LE(w.back(), 0.5);
        } else {
            EXPECT_GT(w.back(), 0.5);
        }
    }
    obs.compute_rows(ranks, q.v.data(), 1);
    // weighted score = sum over quartets of weight * sum_k q log(e / q)
    Network net = nw("((a:1,b:1):1,(c:1,d:1):1,(e:1,f:1):1);").deepcopynetwork();
    semidirectnetwork(net);
    ClassOptions co;
    co.weights = w.data();
    co.perquartet = true;
    co.qlogq = obsqlogq(q, w.data());
    QuartetClasses cls = buildclasses(net, QuartetList::allquartets(6), q, co);
    std::vector<double> params = pidparams(net, cls.npid);
    double want = 0.0;
    for (size_t j = 0; j < static_cast<size_t>(nq); ++j) {
        auto e = classexpectedCF(*cls.eqn[cls.qclass[j]], params.data(), cls.npid, rhotoalpha(0.0));
        for (int k = 0; k < 3; ++k) {
            const double qk = q(j, (cls.qcols[j] >> (2 * k)) & 3);
            if (w[j] > 0 && qk > 0) want += w[j] * qk * std::log(e[static_cast<size_t>(k)] / qk);
        }
    }
    EXPECT_NEAR(scoreclasses(net, cls, 0.0), want, 1e-10);

    // complete gene trees: all weights are 1, and the search is unchanged
    std::vector<Network> gts = read_trees("lazy_gts.tre");
    Network start = read_trees("lazy_tre.tre").at(0);
    double scores[2];
    for (int k = 0; k < 2; ++k) {
        SnaqOptions o = quiet_options();
        o.threads = 1;
        o.hmax = 1;
        o.runs = 2;
        o.Nfail = 10;
        o.geneWeights = k == 1;
        scores[k] = snaq({start}, true, gts, o).fscore;
    }
    EXPECT_EQ(scores[0], scores[1]);
}

TEST(ScoreRank, ClosestToZeroIsBest) {
    // Networks are compared by the distance of their score from 0. Scores are
    // at most 0 up to rounding, where this is the order of the scores.
    EXPECT_EQ(scorerank(-1.5), -1.5);
    EXPECT_EQ(scorerank(0.0), 0.0);
    EXPECT_GT(scorerank(-0.1), scorerank(-0.2));
    EXPECT_GT(scorerank(-0.1), scorerank(0.5));
    EXPECT_GT(scorerank(0.1), scorerank(-0.5));
    EXPECT_EQ(scorerank(0.3), scorerank(-0.3));
}
