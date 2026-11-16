// Ports of SNaQ.jl test/test_inplace_updates/test_rNNI_inplace.jl,
// test/test_perfectData.jl and test/test_negativeEdgeLengths.jl
#include <cmath>

#include "optimize.hpp"
#include "test_support.hpp"

using namespace snaqpp;
using namespace snaqpp::test;

namespace {
void inplace_matches_rebuild(int type) {
    Network net = read_trees("inplace.net").at(0);
    CFMatrix q = read_cf("inplace.cf");
    EquationSet es = findquartetequations(net);
    auto stuvs = allvalidrNNInodes(net, type);
    ASSERT_FALSE(stuvs.empty());
    Rng rng(42);
    for (size_t j = 0; j < stuvs.size(); ++j) {
        // in place
        Network c1 = net.clone();
        auto v1 = allvalidrNNInodes(c1, type)[j];
        performrNNI(c1, v1[0], v1[1], v1[2], v1[3], type, rng);
        OptInfo i1 = gatheroptimizationinfo(c1, true);
        if (!sameparameterlayout(net, c1)) continue;  // search rebuilds the equations of such moves
        // SNaQ.jl's in-place rNNI(1) update misses quartets when v is the
        // parent of u (t's subtree moves below u): SNaQ++'s search rebuilds
        // those, so only the case u = parent of v is compared here.
        if (type == 1 && !(c1.getparents(v1[3]).size() == 1 && c1.getparent_node(v1[3]) == v1[2])) continue;
        std::vector<QuartetData> inplace;
        // SNaQ.jl's test passes the parameter map of the network *before* the
        // move, which raises a KeyError in SNaQ.jl too whenever the move turns
        // a pendant edge into an internal one (the test is not part of
        // SNaQ.jl's runtests.jl). `search` passes the map of the network
        // after the move, as done here.
        if (type == 1)
            applyrNNI1update(c1, es.eqns, inplace, i1.param_map, v1[2], 1);
        else
            applyrNNI2update(c1, es.eqns, inplace, i1.param_map, v1[0], v1[1], v1[2], v1[3], 1);
        double l_in = computeSNaQscore(inplace, i1.params, q, 0.0);
        // rebuilt
        Network c2 = net.clone();
        auto v2 = allvalidrNNInodes(c2, type)[j];
        performrNNI(c2, v2[0], v2[1], v2[2], v2[3], type, rng);
        EquationSet es2 = findquartetequations(c2);
        double l_out = computeSNaQscore(es2.eqns, es2.info.params, q, 0.0);
        EXPECT_EQ(l_in, l_out) << "move " << j;
    }
}
}  // namespace

TEST(InplaceUpdates, rNNI1) { inplace_matches_rebuild(1); }
TEST(InplaceUpdates, rNNI2) { inplace_matches_rebuild(2); }

TEST(PerfectData, ZeroScoreWithExpectedCFs) {
    Network net = nw("(15,(1,((14,(#H1,(((12,13),(11,#H3)),(7,((10)#H3,(8,9)))))),((((2,3))#H2,(6,(5,(#H2,4)))))#H1)));");
    for (EdgeId e : net.edge) net.E(e).length = 0.25;
    for (NodeId h : net.hybrid) {
        net.E(net.getparentedge(h)).gamma = 0.35;
        net.E(net.getparentedgeminor(h)).gamma = 0.65;
    }
    CFMatrix df = computeexpectedCFmatrix(net);
    // optimizetopology!(net, df): equations of all quartets, maxeval = 500
    EquationSet es = findquartetequations(net);
    FitOptions fo;
    fo.maxeval = 500;
    double val = fitnumericalparameters(net, es.eqns, df, 0.0, fo);
    EXPECT_NEAR(val, 0.0, 1e-10);
}

TEST(NegativeEdgeLengths, Inputs) {
    Network net = nw("(((a,#H2),(b,#H1)),((c,d),((e,(f)#H1))#H2));");
    gatheroptimizationinfo(net, true);  // sets the edge numbers so that sett works as expected
    auto sett = [&](EdgeId e) { return net.E(e).hybrid ? 0.0 : std::abs(net.E(e).number / 10.0); };
    for (EdgeId e : net.edge) {
        net.E(e).length = sett(e);
        net.E(e).gamma = net.E(e).hybrid ? 0.5 : -1;
    }
    Network n0 = net.clone();
    CFMatrix dcf0 = computeexpectedCFmatrix(n0);
    for (EdgeId e : net.edge) net.E(e).length = net.N(net.getchild(e)).leaf ? -1 : sett(e);
    Network n1 = net.clone();
    CFMatrix dcf1 = computeexpectedCFmatrix(n1);
    for (EdgeId e : net.edge) {
        NodeId c = net.getchild(e);
        net.E(e).length = (net.N(c).leaf || net.N(c).hybrid) ? -1 : sett(e);
    }
    Network n2 = net.clone();
    CFMatrix dcf2 = computeexpectedCFmatrix(n2);
    ASSERT_EQ(dcf0.v.size(), dcf1.v.size());
    for (size_t i = 0; i < dcf0.v.size(); ++i) {
        EXPECT_EQ(dcf0.v[i], dcf1.v[i]);
        EXPECT_EQ(dcf1.v[i], dcf2.v[i]);
    }
    // a negative length that is not allowed still errors
    net.E(net.edge[9]).length = -1.0;
    Network bad = net.clone();
    EXPECT_THROW(computeexpectedCFmatrix(bad), SnaqError);
    EXPECT_EQ(getnegativeedges(net).size(), 1u);
    // such a network can be used in snaq!
    net.E(net.edge[9]).length = 2.0;
    std::vector<std::string> taxa{"a", "b", "c", "d", "e", "f"};
    SnaqOptions o = quiet_options();
    o.hmax = 2;
    o.runs = 100;
    o.Nfail = 10;
    Network snaqnet = snaq_cfs({net}, true, dcf0, taxa, o);
    EXPECT_NEAR(snaqnet.fscore, 0.0, 1e-6);
}
