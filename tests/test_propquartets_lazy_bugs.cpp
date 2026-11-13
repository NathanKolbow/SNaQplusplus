// Ports of SNaQ.jl test/test_propQuartets.jl, the in-scope parts of
// test/test_lazy_datacf.jl and test/test_bugs.jl (see TESTS_SKIPPED.md).
#include <cmath>
#include <filesystem>
#include <fstream>
#include <regex>
#include <sstream>

#include "quartets.hpp"
#include "search.hpp"
#include "test_support.hpp"

using namespace csnaq;
using namespace csnaq::test;

namespace {

std::string slurp(const std::string& path) {
    std::ifstream in(path);
    std::stringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

struct Data {
    CFMatrix q;
    std::vector<std::string> taxa;
};

Data expected(const Network& net) {
    Data d;
    Network c = net.clone();
    d.q = computeexpectedCFmatrix(c);
    d.taxa = net.tiplabels();
    std::sort(d.taxa.begin(), d.taxa.end());
    return d;
}

double score(const Network& net, const CFMatrix& q) {
    Network c = net.clone();
    return computeSNaQscore_net(c, q);
}

// CFs of all quartets from gene trees, as DataCF(gts) (lazy = false) or
// LazyQuartetCF (lazy = true) computes them
Data from_genetrees(const std::vector<Network>& gts, bool lazy) {
    Data d;
    d.taxa = sorted_union_taxa(gts);
    ObservedCFs obs(gts, d.taxa, lazy);
    std::vector<long long> ranks(static_cast<size_t>(obs.nquartets()));
    for (size_t k = 0; k < ranks.size(); ++k) ranks[k] = static_cast<long long>(k);
    d.q.v.resize(3 * ranks.size());
    obs.compute_rows(ranks, d.q.v.data(), 1);
    return d;
}

const char* kTrue = "(((a,b):0.5,(c,d):0.5):0.5,(e,f):0.5);";
const char* kFalse = "(((a,c):5.0,(e,b):5.0):1.0,(d,f):1.5);";

}  // namespace

// ------------------------------------------------------------ propQuartets

TEST(PropQuartets, SearchImprovesAndFinalScoreUsesAllQuartets) {
    Network truetre = nw(kTrue), falsetre = nw(kFalse);
    Data d = expected(truetre);
    SnaqOptions o = quiet_options();
    o.hmax = 0;
    o.propQuartets = 0.5;
    o.seed = 42;
    o.Nfail = 10;
    o.runs = 100;
    o.maxeval = 5000;
    Network snaqtre = snaq_cfs({falsetre}, true, d.q, d.taxa, o);
    EXPECT_EQ(hardwiredclusterdistance(snaqtre, truetre, false), 0);
    EXPECT_NEAR(snaqtre.fscore, score(truetre, d.q), 1e-10);
}

TEST(PropQuartets, DifferentQuartetsEachRun) {
    Rng r1(1), r2(2);
    EXPECT_NE(sampleqindices(15, 0.5, r1), sampleqindices(15, 0.5, r2));
}

TEST(PropQuartets, FinalZeroSkipsReoptimization) {
    Network truetre = nw(kTrue), falsetre = nw(kFalse);
    Data d = expected(truetre);
    std::string dir = tempdir();
    SnaqOptions o = quiet_options();
    o.hmax = 0;
    o.propQuartets = 0.5;
    o.propQuartetsFinal = 0.0;
    o.seed = 42;
    o.Nfail = 10;
    o.runs = 10;
    o.filename = dir + "/pqf";
    Network snaqtre = snaq_cfs({falsetre}, true, d.q, d.taxa, o);
    EXPECT_EQ(hardwiredclusterdistance(snaqtre, truetre, false), 0);
    EXPECT_EQ(slurp(dir + "/pqf_runs/run1.log").find("Re-optimizing"), std::string::npos);
    std::filesystem::remove_all(dir);
}

TEST(PropQuartets, FinalBetweenZeroAndOneRescoresOnSharedSample) {
    Network truetre = nw(kTrue), falsetre = nw(kFalse);
    Data d = expected(truetre);
    std::string dir = tempdir();
    SnaqOptions o = quiet_options();
    o.propQuartets = 0.5;
    o.propQuartetsFinal = 0.25;
    o.runs = 3;
    o.seed = 42;
    o.Nfail = 10;
    o.filename = dir + "/pqf";
    // multisearch's defaults
    o.ftolRel = o.ftolAbs = o.xtolRel = o.xtolAbs = 1e-12;
    o.probST = 0.3;
    CFSource src = CFSource::full(d.q, d.taxa);
    auto nets = multisearch({falsetre}, true, src, 0, o, 0.25);
    Rng rng(42);
    auto idxs = sampleqindices(nchoose4(6), 0.25, rng);
    CFMatrix Q = src.rows(idxs, 1);
    for (Network& net : nets) {
        Network c = net.clone();
        EquationSet es = findquartetequations(c, idxs);
        EXPECT_NEAR(net.fscore, computeSNaQscore(es.eqns, es.info.params, Q, 0.0), 1e-10);
    }
    std::string run1 = slurp(dir + "/pqf_runs/run1.log"), mainlog = slurp(dir + "/pqf.log");
    EXPECT_NE(run1.find("propQuartetsFinal = 0.25"), std::string::npos);
    EXPECT_NE(mainlog.find("Re-scored the networks of all 3 runs"), std::string::npos);
    EXPECT_EQ(run1.find("Re-scored"), std::string::npos);
    std::filesystem::remove_all(dir);
}

TEST(PropQuartets, FinalOneSameAsDefault) {
    Network truetre = nw(kTrue), falsetre = nw(kFalse);
    Data d = expected(truetre);
    SnaqOptions o = quiet_options();
    o.hmax = 0;
    o.propQuartets = 0.5;
    o.seed = 42;
    o.Nfail = 10;
    o.runs = 5;
    Network a = snaq_cfs({falsetre}, true, d.q, d.taxa, o);
    o.propQuartetsFinal = 1.0;
    Network b = snaq_cfs({falsetre}, true, d.q, d.taxa, o);
    EXPECT_EQ(writenewick(a), writenewick(b));
    EXPECT_EQ(a.fscore, b.fscore);
}

TEST(PropQuartets, FinalMustBeInZeroOne) {
    Network truetre = nw(kTrue), falsetre = nw(kFalse);
    Data d = expected(truetre);
    SnaqOptions o = quiet_options();
    o.hmax = 0;
    o.runs = 1;
    o.propQuartetsFinal = -0.1;
    EXPECT_ANY_THROW(snaq_cfs({falsetre}, true, d.q, d.taxa, o));
    o.propQuartetsFinal = 1.1;
    EXPECT_ANY_THROW(snaq_cfs({falsetre}, true, d.q, d.taxa, o));
}

TEST(PropQuartets, LazyDataCF) {
    Network truetre = nw(kTrue), falsetre = nw(kFalse);
    auto gts = read_trees("propq_gts.tre");
    ASSERT_EQ(gts.size(), 1000u);
    SnaqOptions o = quiet_options();
    o.hmax = 0;
    o.lazy = true;
    o.propQuartets = 1.0;  // a lazy DataCF requires propQuartets < 1
    EXPECT_THROW(snaq({falsetre}, true, gts, o), SnaqError);
    o.propQuartetsFinal = 0.5;
    EXPECT_THROW(snaq({falsetre}, true, gts, o), SnaqError);

    o.propQuartets = 0.8;
    o.propQuartetsFinal = 1.0;
    o.seed = 42;
    o.Nfail = 10;
    o.runs = 10;
    Network snaqtre = snaq({falsetre}, true, gts, o);
    EXPECT_EQ(hardwiredclusterdistance(snaqtre, truetre, false), 0);
    Data lazy = from_genetrees(gts, true), full = from_genetrees(gts, false);
    EXPECT_NEAR(snaqtre.fscore, score(snaqtre, lazy.q), 1e-10);
    EXPECT_NEAR(score(truetre, full.q), score(truetre, lazy.q), 1e-10);
}

// --------------------------------------------------------- lazy DataCF

TEST(LazyDataCF, MatchesNonLazy) {
    Network tre = read_trees("lazy_tre.tre").at(0);
    auto gts = read_trees("lazy_gts.tre");
    Data lazy = from_genetrees(gts, true), full = from_genetrees(gts, false);
    EXPECT_EQ(lazy.q.nrow(), 70u);
    EXPECT_EQ(lazy.taxa, (std::vector<std::string>{"t1", "t2", "t3", "t4", "t5", "t6", "t7", "t8"}));
    EXPECT_NEAR(score(tre, lazy.q), score(tre, full.q), 1e-10);
    ObservedCFs obs(gts, lazy.taxa, true);
    auto q0 = unrank4taxa(8, 0);
    EXPECT_EQ(q0, (std::array<int, 4>{0, 1, 2, 3}));
    EXPECT_THROW(obs.compute(70), SnaqError);
}

// --------------------------------------------------------------- bugs

namespace {
Data bucky() {
    Data d;
    d.q = read_cf("bucky.cf");
    std::istringstream ss(read_lines("bucky_taxa.txt").at(0));
    std::string t;
    while (ss >> t) d.taxa.push_back(t);
    return d;
}

std::string restrictionline(const Restrictions& r) {
    Data d = bucky();
    std::string dir = tempdir();
    SnaqOptions o = quiet_options();
    o.hmax = 0;
    o.runs = 1;
    o.Nfail = 1;
    o.seed = 1;
    o.filename = dir + "/snaq";
    o.restrictions = r;
    snaq_cfs({nw("((((D,C),(A,B)),E),O);")}, true, d.q, d.taxa, o);
    std::string log = slurp(dir + "/snaq.log");
    std::filesystem::remove_all(dir);
    std::smatch m;
    std::regex rx("restrictions = (.*)");
    return std::regex_search(log, m, rx) ? m[1].str() : "";
}
}  // namespace

TEST(Bugs, LogFileDescribesRestrictions) {
    std::string none = restrictionline(Restrictions{});
    EXPECT_EQ(none.find("custom"), std::string::npos);
    EXPECT_NE(none.find("none"), std::string::npos);
    Restrictions g;
    g.galled_tree = true;
    EXPECT_NE(restrictionline(g).find("galled tree"), std::string::npos);
    Restrictions s;
    s.strongly_tree_child = true;
    EXPECT_NE(restrictionline(s).find("tree-child"), std::string::npos);
    Restrictions l;
    l.max_level = 2;
    EXPECT_NE(restrictionline(l).find("level ≤ 2"), std::string::npos);
    Restrictions lg;
    lg.max_level = 2;
    lg.galled_network = true;
    std::string line = restrictionline(lg);
    EXPECT_NE(line.find("level ≤ 2"), std::string::npos);
    EXPECT_NE(line.find("galled network"), std::string::npos);
    Restrictions c;
    c.custom = [](Network&) { return true; };
    c.custom_name = "mycheck";
    EXPECT_NE(restrictionline(c).find("mycheck"), std::string::npos);
}

TEST(Bugs, HmaxEnforcedOnStartingTopologies) {
    Data d = bucky();
    Network starth2 = nw("(((((D,C))#H2,((A,B))#H1),(#H1,E)),(#H2,O));");
    EXPECT_EQ(starth2.numhybrids, 2);
    for (int hmax : {0, 1}) {
        SnaqOptions o = quiet_options();
        o.hmax = hmax;
        o.runs = 2;
        o.Nfail = 5;
        o.seed = 42;
        Network est = snaq_cfs({starth2}, true, d.q, d.taxa, o);
        EXPECT_LE(est.numhybrids, hmax);
    }
    EXPECT_EQ(starth2.numhybrids, 2);
    Network starth2gam =
        nw("(((((D:1,C:1):1)#H2:1::0.8,((A:1,B:1):1)#H1:1::0.7):1,(#H1:1::0.3,E:1):1):1,(#H2:1::0.2,O:1):1);");
    SnaqOptions o = quiet_options();
    CFSource src = CFSource::full(d.q, d.taxa);
    try {
        search(starth2gam, src, 0, o, 1, 1.0, "", "", 1);
        FAIL() << "expected an error";
    } catch (const SnaqError& e) {
        EXPECT_NE(std::string(e.what()).find("hmax"), std::string::npos);
    }
}
