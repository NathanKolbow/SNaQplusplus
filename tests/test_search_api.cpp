// Port of SNaQ.jl test/test_gradient_opt/test_search_API.jl
#include <cmath>
#include <sstream>

#include "compare.hpp"
#include "search.hpp"
#include "snaq.hpp"
#include "test_support.hpp"

using namespace csnaq;
using namespace csnaq::test;

namespace {

struct Truth {
    Network net;
    CFMatrix q;
    std::vector<std::string> taxa;
};

Truth n1_truth() {
    Truth t;
    t.net = read_trees("n1.netfile").at(0);
    Network c = t.net.clone();
    t.q = computeexpectedCFmatrix(c);
    t.taxa = t.net.tiplabels();
    std::sort(t.taxa.begin(), t.taxa.end());
    return t;
}

Truth expected_truth(const Network& net, double rho = 0.0) {
    Truth t;
    t.net = net.clone();
    Network c = net.clone();
    t.q = computeexpectedCFmatrix(c, rho);
    t.taxa = net.tiplabels();
    std::sort(t.taxa.begin(), t.taxa.end());
    return t;
}

double score(const Network& net, const CFMatrix& q, double rho = 0.0) {
    Network c = net.clone();
    return computeSNaQscore_net(c, q, rho);
}

std::vector<Network> multisearch_q(std::vector<Network> starts, const CFMatrix& q, const std::vector<std::string>& taxa,
                                   int hmax, SnaqOptions o) {
    CFSource src = CFSource::full(q, taxa);
    bool single = starts.size() == 1;
    return multisearch(starts, single, src, hmax, o, o.propQuartetsFinal);
}

SnaqOptions msearch_opts(long long seed, int runs) {
    SnaqOptions o = quiet_options();
    o.seed = seed;
    o.runs = runs;
    o.Nfail = 100;  // maxequivPLs = 100
    o.propQuartetsFinal = 1.0;
    // multisearch's defaults differ from snaq!'s for these:
    o.ftolRel = o.ftolAbs = o.xtolRel = o.xtolAbs = 1e-12;
    o.probST = 0.3;
    return o;
}

bool root_has_child_named(const Network& net, const std::string& name) {
    for (NodeId c : net.getchildren(net.getroot()))
        if (net.N(c).name == name) return true;
    return false;
}

}  // namespace

TEST(SearchAPI, MultisearchImprovesOverStart) {
    Truth t = n1_truth();
    Network tre0 = read_trees("search_tre0.tre").at(0);
    auto nets = multisearch_q({tre0.clone()}, t.q, t.taxa, t.net.numhybrids, msearch_opts(7, 10));
    EXPECT_GT(score(nets[0], t.q), score(tre0, t.q));

    auto all = read_trees("search_t0s.tre");
    std::vector<Network> t0s{all[0].clone(), all[1].clone(), all[2].clone()};
    auto nets2 = multisearch_q(t0s, t.q, t.taxa, t.net.numhybrids, msearch_opts(5, 3));
    double minstart = 1e300;
    for (const auto& s : t0s) minstart = std::min(minstart, score(s, t.q));
    EXPECT_GT(score(nets2[0], t.q), minstart);
}

TEST(SearchAPI, Outgroups) {
    Truth t = n1_truth();
    auto all = read_trees("search_t0s.tre");
    const char* outgroups[] = {"a", "b", "e", "f"};
    for (int k = 0; k < 4; ++k) {
        std::vector<Network> t0s{all[3 + 3 * k].clone(), all[4 + 3 * k].clone(), all[5 + 3 * k].clone()};
        SnaqOptions o = msearch_opts(5, 3);
        o.outgroup = outgroups[k];
        auto nets = multisearch_q(t0s, t.q, t.taxa, t.net.numhybrids, o);
        EXPECT_TRUE(root_has_child_named(nets[0], outgroups[k])) << outgroups[k];
        if (k >= 2) EXPECT_NE(hardwiredclusterdistance(t.net, nets[0], false), 0) << outgroups[k];
    }
}

TEST(SearchAPI, SnaqWithPropQuartets) {
    Truth t = n1_truth();
    auto t0s = read_trees("search_t0.tre");
    {
        Network t0 = t0s[0].clone();
        double initLL = score(t0, t.q);
        SnaqOptions o = quiet_options();
        o.hmax = 1;
        o.propQuartets = 0.85;
        o.runs = 5;
        o.Nfail = 50;
        Network snaqnet = snaq_cfs({t0}, true, t.q, t.taxa, o);
        EXPECT_NEAR(snaqnet.fscore, score(snaqnet, t.q), 1e-8);
        EXPECT_GT(snaqnet.fscore, initLL);
    }
    {
        Network t0 = t0s[1].clone();
        double initLL = score(t0, t.q);
        SnaqOptions o = quiet_options();
        o.hmax = 1;
        o.propQuartets = 0.85;
        o.runs = 5;
        o.Nfail = 50;
        o.restrictions.max_level = 1;
        Network snaqnet = snaq_cfs({t0}, true, t.q, t.taxa, o);
        EXPECT_TRUE(o.restrictions(snaqnet));
        EXPECT_TRUE(tcgidentifiable(snaqnet));
        EXPECT_GT(snaqnet.fscore, initLL);
        EXPECT_NEAR(snaqnet.fscore, score(snaqnet, t.q), 1e-8);
    }
}

TEST(SearchAPI, CustomRestrictions) {
    Truth t = n1_truth();
    Network t0 = read_trees("search_t0.tre").at(2);
    SnaqOptions o = quiet_options();
    o.restrictions.custom = [](Network& n) { return tcgidentifiable(n); };
    o.restrictions.custom_name = "tcgidentifiable";
    o.hmax = 4;
    o.runs = 20;
    o.Nfail = 50;
    Network snaqnet = snaq_cfs({t0}, true, t.q, t.taxa, o);
    EXPECT_TRUE(tcgidentifiable(snaqnet));
}

TEST(SearchAPI, MaxLevelZeroInfersTree) {
    Truth t = n1_truth();
    Network t0 = read_trees("search_t0.tre").at(3);
    SnaqOptions o = quiet_options();
    o.restrictions.max_level = 0;
    o.hmax = 4;
    o.runs = 20;
    o.Nfail = 50;
    Network snaqnet = snaq_cfs({t0}, true, t.q, t.taxa, o);
    EXPECT_EQ(snaqnet.numhybrids, 0);
    EXPECT_TRUE(tcgidentifiable(snaqnet));
}

TEST(SearchAPI, RespectsRestrictionsWhenTruthOutside) {
    Network truenet = nw("(((t4:0.4,(t8:0.3,t2:0.3):0.1):0.2,(t1:0.4,((t5:0.2,#H2:0.0::0.25):0.2,((t6:0.0,(t3:0.0)#H1:0.0::0.75):0.2)#H2:0.2::0.75):0.0):0.2):0.6,(t7:0.6,#H1:0.0::0.25):0.6);");
    Truth t = expected_truth(truenet);
    Network T0 = read_trees("search_level2_T0.tre").at(0);
    SnaqOptions o = quiet_options();
    o.restrictions.max_level = 1;
    o.hmax = 2;
    o.runs = 2;
    o.Nfail = 25;
    Network snaqnet = snaq_cfs({T0.clone()}, true, t.q, t.taxa, o);
    EXPECT_GT(score(snaqnet, t.q), score(T0, t.q));
    EXPECT_EQ(getnetworklevel(snaqnet), 1);
}

TEST(SearchAPI, ManyHybridsDoesNotErrorAndImproves) {
    auto lines = read_lines("search_manyhybrids.txt");
    ASSERT_EQ(lines.size(), 180u);
    for (size_t k = 0; k < lines.size(); k += 3) {
        std::istringstream ss(lines[k]);
        int n, h, seed;
        ss >> n >> h >> seed;
        Network truenet = nw(lines[k + 1]);
        Network T0 = nw(lines[k + 2]);
        Truth t = expected_truth(truenet);
        SnaqOptions o = quiet_options();
        o.restrictions.custom = [](Network& x) { return tcgidentifiable(x); };
        o.restrictions.custom_name = "tcgidentifiable";
        o.propQuartets = 0.1;
        o.hmax = h;
        o.runs = 2;
        o.Nfail = 5;
        Network snaqnet;
        ASSERT_NO_THROW(snaqnet = snaq_cfs({T0.clone()}, true, t.q, t.taxa, o)) << lines[k];
        EXPECT_GT(score(snaqnet, t.q), score(T0, t.q)) << lines[k];
        EXPECT_TRUE(tcgidentifiable(snaqnet)) << lines[k];
    }
}

TEST(SearchAPI, MissingEdgeLengthsAndGammasInStart) {
    Network net = read_trees("gennet_10_2_8.net").at(0);
    Truth t = expected_truth(net);
    for (EdgeId e : net.edge) net.E(e).length = -1;
    for (NodeId h : net.hybrid) {
        net.E(net.getparentedge(h)).gamma = -1;
        net.E(net.getparentedgeminor(h)).gamma = -1;
    }
    SnaqOptions o = quiet_options();
    o.hmax = 2;
    o.Nfail = 20;
    o.runs = 2;
    o.opt_maxeval = 1000;
    o.probST = 0.0;
    Network snaqnet = snaq_cfs({net}, true, t.q, t.taxa, o);
    EXPECT_GT(score(snaqnet, t.q), -1e-2);
}

TEST(SearchAPI, PolytomiesInInput) {
    {
        Network multi = nw("((a,b),(c,d),(e,f),(g,(h,i,j)));");
        for (EdgeId e : multi.edge) multi.E(e).length = 1.0;
        Truth t = expected_truth(multi);
        for (EdgeId e : multi.edge) multi.E(e).length = -1;
        SnaqOptions o = quiet_options();
        o.hmax = 0;
        o.Nfail = 20;
        o.runs = 2;
        o.probST = 0.75;
        Network snaqnet = snaq_cfs({multi.clone()}, true, t.q, t.taxa, o);
        std::vector<EdgeId> rev(snaqnet.edge.rbegin(), snaqnet.edge.rend());
        for (EdgeId e : rev) {
            if (snaqnet.index_of_edge(e) < 0) continue;
            if (!snaqnet.N(snaqnet.getchild(e)).leaf && snaqnet.E(e).length >= 0.0 && snaqnet.E(e).length <= 1e-5)
                shrinkedge(snaqnet, e);
        }
        EXPECT_EQ(hardwiredclusterdistance(snaqnet, multi, false), 0);
    }
    {
        Network multi = read_trees("multifurcation_net.net").at(0);
        Truth t = expected_truth(multi);
        SnaqOptions o = quiet_options();
        o.hmax = multi.numhybrids;
        o.Nfail = 20;
        o.runs = 2;
        o.opt_maxeval = 1000;
        Network snaqnet = snaq_cfs({multi.clone()}, true, t.q, t.taxa, o);
        EXPECT_GT(snaqnet.fscore, -1e-6);
        std::vector<EdgeId> rev(snaqnet.edge.rbegin(), snaqnet.edge.rend());
        for (EdgeId e : rev) {
            if (snaqnet.index_of_edge(e) < 0) continue;
            if (!snaqnet.N(snaqnet.getchild(e)).leaf && !snaqnet.E(e).hybrid && snaqnet.E(e).length >= 0.0 &&
                snaqnet.E(e).length <= 1e-5)
                shrinkedge(snaqnet, e);
        }
        EXPECT_EQ(hardwiredclusterdistance(snaqnet, multi, false), 0);
    }
}

namespace {
void start_vs_search(SnaqOptions o, double rho_data, double rho_score) {
    Network tnet = read_trees("gennet_10_2_8.net").at(0);
    Truth t = expected_truth(tnet, rho_data);
    Network startnet = read_trees("gennet_10_2_12.net").at(0);
    double startLL = score(startnet, t.q, rho_score);
    o.hmax = tnet.numhybrids;
    o.Nfail = 20;
    o.runs = 2;
    o.opt_maxeval = 100;
    Network snaqnet = snaq_cfs({startnet.clone()}, true, t.q, t.taxa, o);
    EXPECT_GE(score(snaqnet, t.q, rho_score), startLL);
}
}  // namespace

TEST(SearchAPI, ProbQRHalf) {
    SnaqOptions o = quiet_options();
    o.probQR = 0.5;
    o.probST = 0.0;
    start_vs_search(o, 0.0, 0.0);
}

TEST(SearchAPI, ProbQROne) {
    SnaqOptions o = quiet_options();
    o.probQR = 1.0;
    o.probST = 0.0;
    start_vs_search(o, 0.0, 0.0);
}

TEST(SearchAPI, QinfTest) {
    SnaqOptions o = quiet_options();
    o.qinfTest = true;
    o.probST = 0.0;
    start_vs_search(o, 0.0, 0.0);
}

TEST(SearchAPI, QinfTestProbQRProbST) {
    SnaqOptions o = quiet_options();
    o.qinfTest = true;
    o.probQR = 0.9;
    o.probST = 1.0;
    start_vs_search(o, 0.0, 0.0);
}

TEST(SearchAPI, QinfTestProbQRPropQuartetsProbST) {
    SnaqOptions o = quiet_options();
    o.qinfTest = true;
    o.propQuartets = 0.25;
    o.probQR = 0.9;
    o.probST = 1.0;
    start_vs_search(o, 0.0, 0.0);
}

TEST(SearchAPI, RhoOne) {
    SnaqOptions o = quiet_options();
    o.rho = 1.0;
    o.probST = 0.0;
    start_vs_search(o, 0.0, 0.0);
}

TEST(SearchAPI, RhoOneDataRhoOne) {
    SnaqOptions o = quiet_options();
    o.rho = 1.0;
    o.probST = 0.0;
    start_vs_search(o, 1.0, 1.0);
}

// snaq --optimize: only the parameters of the given topology change.
TEST(OptimizeNetwork, FitsParametersOfFixedTopology) {
    std::vector<Network> gts = read_trees("lazy_gts.tre");
    Network start = read_trees("lazy_tre.tre").at(0);
    SnaqOptions o;
    o.seed = 3;
    o.seed_given = true;
    Network fitted = optimize_network(start.clone(), gts, o);
    EXPECT_TRUE(std::isfinite(fitted.fscore));
    EXPECT_LE(fitted.fscore, 1e-6);
    EXPECT_EQ(hardwiredclusterdistance(fitted, start, false), 0);

    // a sample of the quartets
    o.propQuartets = 0.3;
    Network sub = optimize_network(start.clone(), gts, o);
    EXPECT_TRUE(std::isfinite(sub.fscore));
    EXPECT_EQ(hardwiredclusterdistance(sub, start, false), 0);
    EXPECT_THROW(
        {
            SnaqOptions bad = o;
            bad.propQuartets = 0.0;
            optimize_network(start.clone(), gts, bad);
        },
        SnaqError);
}
