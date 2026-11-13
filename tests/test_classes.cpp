// Tests of the quartet-class likelihood (classes.hpp): expected CFs against
// an independent brute-force calculation, gradients against finite
// differences, and equations reused after topology moves against equations
// built from scratch. (Not ported from SNaQ.jl.)
#include <cmath>

#include "bruteforce_cf.hpp"
#include "classes.hpp"
#include "porteqn.hpp"
#include "test_support.hpp"

using namespace csnaq;
using namespace csnaq::test;

namespace {

void randomize_parameters(Network& net, Rng& rng) {
    for (EdgeId e : net.edge) net.E(e).length = 0.05 + 1.5 * rng.rand();
    for (NodeId h : net.hybrid) {
        const double g = 0.05 + 0.4 * rng.rand();
        net.E(net.getparentedgeminor(h)).gamma = g;
        net.E(net.getparentedge(h)).gamma = 1.0 - g;
    }
}

CFMatrix random_cfs(size_t nq, Rng& rng) {
    CFMatrix q;
    q.v.resize(3 * nq);
    for (size_t j = 0; j < nq; ++j) {
        double a = 0.05 + rng.rand(), b = 0.05 + rng.rand(), c = 0.05 + rng.rand();
        q(j, 0) = a / (a + b + c);
        q(j, 1) = b / (a + b + c);
        q(j, 2) = c / (a + b + c);
    }
    return q;
}

std::vector<NodeId> sorted_leaves(const Network& net) {
    std::vector<std::pair<std::string, NodeId>> lv;
    for (NodeId l : net.leaf) lv.emplace_back(net.N(l).name, l);
    std::sort(lv.begin(), lv.end());
    std::vector<NodeId> out;
    for (auto& p : lv) out.push_back(p.second);
    return out;
}

// Expected CFs of every quartet from its class, in the quartet's own column order.
std::vector<std::array<double, 3>> quartet_ecfs(const Network& net, const QuartetClasses& cls, double rho) {
    std::vector<double> params = pidparams(net, cls.npid);
    std::vector<std::array<double, 3>> out(cls.qclass.size());
    for (size_t j = 0; j < out.size(); ++j) {
        auto e = classexpectedCF(*cls.eqn[cls.qclass[j]], params.data(), cls.npid, rhotoalpha(rho));
        for (int k = 0; k < 3; ++k) out[j][(cls.qcols[j] >> (2 * k)) & 3] = e[static_cast<size_t>(k)];
    }
    return out;
}

// Number of quartets whose class expected CFs (ρ = 0) differ from brute force.
long wrong_vs_bruteforce(const Network& net, const QuartetClasses& cls, long* checked = nullptr) {
    BruteForceCF bf(net);
    std::vector<NodeId> leaves = sorted_leaves(net);
    auto ecf = quartet_ecfs(net, cls, 0.0);
    long bad = 0, n = 0;
    for (size_t j = 0; j < ecf.size(); ++j) {
        auto t = unrank4taxa(net.numtaxa, static_cast<long long>(j));
        double want[3];
        if (!bf.ecf({leaves[static_cast<size_t>(t[0])], leaves[static_cast<size_t>(t[1])],
                     leaves[static_cast<size_t>(t[2])], leaves[static_cast<size_t>(t[3])]},
                    want))
            continue;
        ++n;
        if (std::abs(ecf[j][0] - want[0]) + std::abs(ecf[j][1] - want[1]) + std::abs(ecf[j][2] - want[2]) > 1e-9) ++bad;
    }
    if (checked) *checked = n;
    return bad;
}

std::vector<Network> fixture_networks() {
    std::vector<Network> nets;
    nets.push_back(nw("(15,(1,((14,(#H1,(((12,13),(11,#H3)),(7,((10)#H3,(8,9)))))),((((2,3))#H2,(6,(5,(#H2,4)))))#H1)));"));
    for (const char* f : {"gennet_10_2_8.net", "gennet_10_2_12.net"})
        for (Network& n : read_trees(f)) nets.push_back(std::move(n));
    nets.push_back(nw("((a,b),(c,d),((e,f),(g,h)));"));
    for (Network& n : nets) {
        n = n.deepcopynetwork();
        semidirectnetwork(n);
    }
    return nets;
}

QuartetClasses allclasses(Network& net, const CFMatrix& q, const QuartetClasses* reuse = nullptr) {
    ClassOptions co;
    co.perquartet = true;
    co.qlogq = obsqlogq(q);
    co.reuse = reuse;
    return buildclasses(net, QuartetList::allquartets(net.numtaxa), q, co);
}

}  // namespace

TEST(Classes, ExpectedCFsMatchBruteForce) {
    Rng rng(7);
    long total = 0;
    for (Network& net : fixture_networks()) {
        randomize_parameters(net, rng);
        CFMatrix q = random_cfs(static_cast<size_t>(nchoose4(net.numtaxa)), rng);
        QuartetClasses cls = allclasses(net, q);
        EXPECT_LT(cls.size(), q.nrow());
        long checked = 0;
        EXPECT_EQ(wrong_vs_bruteforce(net, cls, &checked), 0);
        total += checked;
    }
    EXPECT_GT(total, 1000);
}

TEST(Classes, ScoreMatchesPerQuartetSum) {
    // The score from class sums equals the sum over quartets of q log(e / q),
    // with each quartet's expected CFs taken from its class (ρ = 0 and 0.5).
    Rng rng(11);
    for (Network& net : fixture_networks()) {
        randomize_parameters(net, rng);
        CFMatrix q = random_cfs(static_cast<size_t>(nchoose4(net.numtaxa)), rng);
        QuartetClasses cls = allclasses(net, q);
        for (double rho : {0.0, 0.5}) {
            auto ecf = quartet_ecfs(net, cls, rho);
            double want = 0.0;
            for (size_t j = 0; j < ecf.size(); ++j)
                for (int k = 0; k < 3; ++k) want += q(j, k) * std::log(ecf[j][static_cast<size_t>(k)] / q(j, k));
            double got = scoreclasses(net, cls, rho);
            EXPECT_NEAR(got, want, 1e-9 * std::max(1.0, std::abs(want)));
        }
    }
}

TEST(Classes, GradientMatchesFiniteDifferences) {
    Rng rng(13);
    for (Network& net : fixture_networks()) {
        randomize_parameters(net, rng);
        CFMatrix q = random_cfs(static_cast<size_t>(nchoose4(net.numtaxa)), rng);
        QuartetClasses cls = allclasses(net, q);
        OptInfo info = gatheroptimizationinfo(net, false);
        for (double rho : {0.0, 0.5}) {
            std::vector<double> params = pidparams(net, cls.npid), grad(params.size()), tmp(params.size());
            LossWorkspace ws;
            evalclasses(cls.batch, params.data(), grad.data(), rhotoalpha(rho), ws, 1);
            for (const ParamObj& o : info.idx_obj) {
                const size_t p = static_cast<size_t>(o.is_hybrid_node ? net.N(o.id).pid : net.E(o.id).pid);
                const double h = 1e-6, x0 = params[p];
                params[p] = x0 + h;
                double fp = evalclasses(cls.batch, params.data(), tmp.data(), rhotoalpha(rho), ws, 1);
                params[p] = x0 - h;
                double fm = evalclasses(cls.batch, params.data(), tmp.data(), rhotoalpha(rho), ws, 1);
                params[p] = x0;
                EXPECT_NEAR(grad[p], (fp - fm) / (2 * h), 1e-5 * std::max(1.0, std::abs(grad[p])));
            }
        }
    }
}

TEST(Classes, MatchPerQuartetEquationsWithInheritanceCorrelation) {
    // Where the per-quartet equations give the same score at ρ = 0, they must
    // also give the same score at ρ = 0.5.
    Rng rng(17);
    int compared = 0;
    for (Network& net : fixture_networks()) {
        randomize_parameters(net, rng);
        CFMatrix q = random_cfs(static_cast<size_t>(nchoose4(net.numtaxa)), rng);
        Network a = net.deepcopynetwork();
        EquationSet es = findquartetequations(a);
        QuartetClasses cls = allclasses(net, q);
        double old0 = computeSNaQscore(es.eqns, es.info.params, q, 0.0), new0 = scoreclasses(net, cls, 0.0);
        if (std::abs(old0 - new0) > 1e-9 * std::abs(old0)) continue;
        ++compared;
        double old5 = computeSNaQscore(es.eqns, es.info.params, q, 0.5), new5 = scoreclasses(net, cls, 0.5);
        EXPECT_NEAR(old5, new5, 1e-9 * std::abs(old5));
    }
    EXPECT_GE(compared, 2);
}

TEST(Classes, EquationsReusedAfterMovesMatchRebuild) {
    Rng rng(23);
    Network N = random_tree(11, rng).deepcopynetwork();
    semidirectnetwork(N);
    N.compact_pids();
    randomize_parameters(N, rng);
    const int hmax = 3;
    CFMatrix q = random_cfs(static_cast<size_t>(nchoose4(N.numtaxa)), rng);
    QuartetClasses N_cls = allclasses(N, q);
    std::vector<AttemptedMove> none;
    std::vector<double> deltas;
    CFMatrix emptyQ;
    int compared = 0, reusedsome = 0, maxh = 0;
    for (int it = 0; it < 400; ++it) {
        Network P = N.deepcopynetwork();
        try {
            auto mp = generatemoveproposal(P, nullptr, none, hmax, 0.0, emptyQ, deltas, rng, 0.0);
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
        randomize_parameters(P, rng);
        QuartetClasses reused = allclasses(P, q, &N_cls);
        Network P2 = P.deepcopynetwork();
        QuartetClasses fresh = allclasses(P2, q);
        ASSERT_EQ(fresh.nbuilt, fresh.size());
        if (reused.nbuilt < reused.size()) ++reusedsome;
        for (double rho : {0.0, 0.5}) {
            auto e1 = quartet_ecfs(P, reused, rho), e2 = quartet_ecfs(P2, fresh, rho);
            ASSERT_EQ(e1.size(), e2.size());
            double worst = 0.0;
            for (size_t j = 0; j < e1.size(); ++j)
                for (size_t k = 0; k < 3; ++k) worst = std::max(worst, std::abs(e1[j][k] - e2[j][k]));
            ASSERT_LT(worst, 1e-12) << "iteration " << it;
        }
        if (it % 10 == 0) ASSERT_EQ(wrong_vs_bruteforce(P, reused), 0) << "iteration " << it;
        ++compared;
        maxh = std::max(maxh, P.numhybrids);
        if (rng.rand() < 0.5) {
            N = std::move(P);
            N_cls = std::move(reused);
            // as after a fit: exchange the major and minor edges of some hybrids
            bool flipped = false;
            for (NodeId h : N.hybrid)
                if (rng.rand() < 0.3) {
                    EdgeId ma = N.getparentedge(h), mi = N.getparentedgeminor(h);
                    N.E(ma).ismajor = false;
                    N.E(mi).ismajor = true;
                    std::swap(N.E(ma).gamma, N.E(mi).gamma);
                    N.E(mi).gamma = 1.0 - N.E(ma).gamma;
                    flipped = true;
                }
            if (flipped) N_cls = allclasses(N, q, &N_cls);
        }
    }
    EXPECT_GT(compared, 200);
    EXPECT_GT(reusedsome, 100);
    EXPECT_GE(maxh, 2);
}

TEST(Classes, EquationBuilderMatchesReference) {
    // PortEquationBuilder (the lineages move up an unchanged network) against
    // the reference implementation (SNaQ.jl's recursion on copies of the
    // network): same expected CFs for every class, and never a larger equation.
    Rng rng(29);
    size_t reticulate = 0;
    auto compare = [&](Network& net, const CFMatrix& q, const std::string& what) {
        ClassOptions co;
        co.qlogq = obsqlogq(q);
        Network a = net.deepcopynetwork(), b = net.deepcopynetwork();
        QuartetClasses fast = buildclasses(a, QuartetList::allquartets(a.numtaxa), q, co);
        co.referenceequations = true;
        QuartetClasses ref = buildclasses(b, QuartetList::allquartets(b.numtaxa), q, co);
        ASSERT_EQ(fast.size(), ref.size()) << what;
        std::vector<double> params = pidparams(a, fast.npid);
        for (size_t ci = 0; ci < fast.size(); ++ci) {
            if (ref.eqn[ci]->eqn.division_H < 0) continue;
            ++reticulate;
            EXPECT_LE(equationsize(fast.eqn[ci]->eqn), equationsize(ref.eqn[ci]->eqn)) << what << ", class " << ci;
            for (double rho : {0.0, 0.5}) {
                auto e1 = classexpectedCF(*fast.eqn[ci], params.data(), fast.npid, rhotoalpha(rho));
                auto e2 = classexpectedCF(*ref.eqn[ci], params.data(), ref.npid, rhotoalpha(rho));
                for (size_t k = 0; k < 3; ++k) ASSERT_NEAR(e1[k], e2[k], 1e-12) << what << ", class " << ci;
            }
        }
    };
    int nfix = 0;
    for (Network& net : fixture_networks()) {
        randomize_parameters(net, rng);
        CFMatrix q = random_cfs(static_cast<size_t>(nchoose4(net.numtaxa)), rng);
        compare(net, q, "fixture " + std::to_string(nfix++));
    }
    // networks reached by random moves, with up to 4 reticulations
    Network N = random_tree(10, rng).deepcopynetwork();
    semidirectnetwork(N);
    N.compact_pids();
    CFMatrix q = random_cfs(static_cast<size_t>(nchoose4(N.numtaxa)), rng);
    std::vector<AttemptedMove> none;
    std::vector<double> deltas;
    CFMatrix emptyQ;
    int compared = 0, maxh = 0;
    for (int it = 0; it < 300; ++it) {
        Network P = N.deepcopynetwork();
        try {
            auto mp = generatemoveproposal(P, nullptr, none, 4, 0.0, emptyQ, deltas, rng, 0.0);
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
        randomize_parameters(P, rng);
        compare(P, q, "move " + std::to_string(it));
        if (::testing::Test::HasFatalFailure()) return;
        ++compared;
        maxh = std::max(maxh, P.numhybrids);
        N = std::move(P);
    }
    EXPECT_GT(compared, 150);
    EXPECT_GE(maxh, 3);
    EXPECT_GT(reticulate, 5000u);
}
