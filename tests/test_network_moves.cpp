// Ports of SNaQ.jl test/test_network_moves/*.jl
#include <cmath>

#include "test_support.hpp"

using namespace csnaq;
using namespace csnaq::test;

namespace {

Network reload_labelled_net(bool named_root = false) {
    Network net = nw(named_root ? "((a,b)i1,((c,#H1)i6,(d,((e,f)i5)#H1)i3)i2)r;"
                                : "((a,b)i1,((c,#H1)i6,(d,((e,f)i5)#H1)i3)i2);");
    semidirectnetwork(net);
    net.N(net.hybrid[0]).name = "i4";
    return net;
}

bool all_gammas_sum_to_one(const Network& net) {
    for (NodeId h : net.hybrid) {
        double s = net.E(net.getparentedge(h)).gamma + net.E(net.getparentedgeminor(h)).gamma;
        if (std::abs(s - 1.0) > 1e-8) return false;
    }
    return true;
}

EdgeId edge_between(const Network& net, NodeId a, NodeId b) {
    for (EdgeId e : net.N(a).edge)
        if (net.E(e).has_node(b)) return e;
    throw SnaqError("not adjacent");
}

std::string replace_all(std::string s, const std::string& from, const std::string& to) {
    size_t p = 0;
    while ((p = s.find(from, p)) != std::string::npos) {
        s.replace(p, from.size(), to);
        p += to.size();
    }
    return s;
}

}  // namespace

// ---------------------------------------------------------------- add/remove

TEST(AddRemoveRetic, InvalidMoves) {
    Network net = nw("(a,((b,#H1)i2,(c,((d,e)i5)#H1)i4)i3)i1;");
    EdgeId e2 = net.getparentedge(net.leaf[3]);
    EdgeId e1 = net.getparentedge(net.getparent(e2));
    EXPECT_TRUE(isvalidaddhybrid(net, e1, e2));
    EXPECT_FALSE(isvalidaddhybrid(net, e2, e1));
}

TEST(AddRemoveRetic, AddIsReversible) {
    Network net = nw("(a,((b,#H1)i2,(c,((e,d)i5)#H1)i4)i3)i1;");
    Network net0 = net.clone();
    EdgeId e1 = net.getparentedge(net.leaf[3]);
    EdgeId e2 = net.getparentedge(net.getparent(e1));
    Rng rng(1);
    NodeId newH = addhybrid(net, e2, e1, rng);
    removehybrid(net, newH);
    EXPECT_EQ(replace_all(writenewick(net), "d", "e"), replace_all(writenewick(net0), "d", "e"));
}

TEST(AddRemoveRetic, MoveOrigin) {
    Network net = nw("(((d1:1.602,d2:1.602)i10:0.823,(b:0.463,(a:0.573)#H1:0.5::0.973)i3442:1.191)i13:0.111,((c2:0.757,c1:1.554)i6:0.799,(e:0.776,f:0.776)i5:0.497)i12:0.0,#H1:0.5::0.027)i7;");
    net.isrooted = false;
    Rng rng(0);
    EXPECT_TRUE(samplemovereticulateoriginparameters(net, rng).valid);
    // the move SNaQ.jl's test samples with Random.seed!(0): H1 to edge i10 -> d2
    auto n = get_nodes(net, {"H1", "i10", "d2"});
    EXPECT_TRUE(isvalidmovereticulateorigin(net, n[0], edge_between(net, n[1], n[2])));
    EXPECT_NO_THROW(movereticulateorigin(net, n[0], edge_between(net, n[1], n[2])));
    EXPECT_EQ(writenewick(net), "((c2:0.757,c1:1.554)i6:0.799,(e:0.776,f:0.776)i5:0.497,((d1:1.602,(d2:0.801,#H1:0.5::0.027)i7:0.801)i10:0.823,(b:0.463,(a:0.573)#H1:0.5::0.973)i3442:1.191)i13:0.111)i12;");
    EdgeId new_origin = net.getparentedge(net.leaf[2]);
    EXPECT_TRUE(isvalidmovereticulateorigin(net, net.hybrid[0], new_origin));

    Network net2 = nw("(b:1.095,(a:0.268,#H2:0.562::0.5)i3649:0.5,((((e:0.679,f:0.515)i5:0.496,(c2:0.51,c1:0.0)i13:0.8)i6:0.109,(d1:0.65,d2:1.602)H1:0.825)i10737:0.562)#H2:0.562::0.5)i7;");
    net2.isrooted = false;
    Rng rng2(0);
    EXPECT_TRUE(samplemovereticulateoriginparameters(net2, rng2).valid);
}

TEST(AddRemoveRetic, RemoveExcessReticulations) {
    Network neth2 = nw("(((((D:1,C:1):1)#H2:1::0.8,((A:1,B:1):1)#H1:1::0.7):1,(#H1:1::0.3,E:1):1):1,(#H2:1::0.2,O:1):1);");
    semidirectnetwork(neth2);
    {
        Network net = neth2.clone();
        removeexcessreticulations(net, 2);
        EXPECT_EQ(net.numhybrids, 2);
        Network a = net.clone(), b = neth2.clone();
        EXPECT_EQ(writenewick(a), writenewick(b));
    }
    {
        Network net = neth2.clone();
        removeexcessreticulations(net, 1);
        ASSERT_EQ(net.numhybrids, 1);
        NodeId H = net.hybrid[0];
        // hardwired cluster below the major parent edge of the kept hybrid
        std::vector<std::string> below;
        std::vector<NodeId> stack{H};
        while (!stack.empty()) {
            NodeId n = stack.back();
            stack.pop_back();
            if (net.N(n).leaf) below.push_back(net.N(n).name);
            for (NodeId c : net.getchildren(n)) stack.push_back(c);
        }
        std::sort(below.begin(), below.end());
        EXPECT_EQ(below, (std::vector<std::string>{"A", "B"}));
        EXPECT_NEAR(net.E(net.getparentedgeminor(H)).gamma, 0.3, 1e-12);
        EXPECT_NEAR(net.E(net.getparentedge(H)).gamma, 0.7, 1e-12);
    }
    {
        Network net = neth2.clone();
        removeexcessreticulations(net, 0);
        EXPECT_EQ(net.numhybrids, 0);
        Network major = read_trees("neth2_majortree.tre").at(0);
        EXPECT_EQ(hardwiredclusterdistance(net, major, false), 0);
    }
}

// ----------------------------------------------------------- origin/target

TEST(MoveTargetOrigin, PastBugs) {
    {
        Network net = nw("(((c1:0.5,c2:0.5)i3:0.622,((e:0.5,f:0.5)i10:0.206)#Hi4:0.206::0.739)i4:0.31,((d2:0.25,(((#Hi4:0.206::0.261,d1:0.5)i7:0.008)#H2:0.008::0.978)#H3:0.008::0.774)i5440:0.662,#H3:0.008::0.226)i1976:0.662,(b:0.5,(a:0.25,#H2:0.008::0.022)i12:0.709)i15:1.004)i16;");
        net.isrooted = false;
        auto n = get_nodes(net, {"H3", "i7", "Hi4"});
        movereticulatetarget(net, n[0], edge_between(net, n[1], n[2]));
        EXPECT_TRUE(all_gammas_sum_to_one(net));
    }
    {
        Network net = nw("(((b:0.75)#H3:0.5::0.972,a:0.96)i15:1.081,(((e:0.5,f:0.5)i10:0.208)#Hi4:0.208::0.72,((d1:0.5,d2:0.75)i7:0.095,#H4:0.5::0.007)i1817:0.549)i12:0.281,((#Hi4:0.208::0.28,(c2:0.447)#H4:0.5::0.993)i3:0.0,(c1:0.132,#H3:0.5::0.028)i5645:0.993)i4:1.127)i16;");
        net.isrooted = false;
        auto n = get_nodes(net, {"Hi4", "i12", "i1817"});
        movereticulateorigin(net, n[0], edge_between(net, n[1], n[2]));
        EXPECT_TRUE(all_gammas_sum_to_one(net));
    }
    {
        Network net = nw("((c2:1.162,(((a:1.136,(b:0.548,#H5:0.5::0.043)i8475:0.625)H3:1.187,((d1:1.038,d2:1.301)i10:0.35,#H4:0.5::0.32)i9744:0.539)i3628:0.149,((f:0.776,(#H6:0.5::0.465)#H4:0.5::0.68)i5:0.481,(e:0.15)#H6:0.5::0.535)i3056:0.481)i12:0.831)i2450:0.5,(c1:0.285)#H5:0.5::0.957)i7;");
        net.isrooted = false;
        auto n = get_nodes(net, {"H4", "i12", "i3056"});
        movereticulatetarget(net, n[0], edge_between(net, n[1], n[2]));
        EXPECT_TRUE(all_gammas_sum_to_one(net));
    }
}

// ------------------------------------------------------------------- rSPR

TEST(rSPR, VerifiedByHand) {
    Network net = reload_labelled_net(true);
    auto v = get_nodes(net, {"i1", "i3", "i6", "i2", "i4", "i5"});
    EXPECT_TRUE(isvalidrSPR(net, v[0], v[1], v[2], v[3], v[4], v[5]));
    auto w = get_nodes(net, {"i1", "i3", "i6", "i2", "i5", "i4"});
    EXPECT_FALSE(isvalidrSPR(net, w[0], w[1], w[2], w[3], w[4], w[5]));
    auto x = get_nodes(net, {"i3", "i1", "i2", "i6", "i4", "i5"});
    EXPECT_FALSE(isvalidrSPR(net, x[0], x[1], x[2], x[3], x[4], x[5]));

    Rng rng(0);
    MoveParams mp;
    mp.valid = true;
    for (NodeId n : v) mp.p.push_back({false, n});
    applymove(net, Move::rSPR, mp, rng);
    EXPECT_TRUE(testnewickequality(net, "(d,(((e,f)i5,(a,b)i1)i2)#i4,(c,#i4)i6)i3;"));

    auto y = get_nodes(net, {"i6", "i3", "i2", "i4", "i5", "e"});
    MoveParams mp2;
    mp2.valid = true;
    for (NodeId n : y) mp2.p.push_back({false, n});
    applymove(net, Move::rSPR, mp2, rng);
    EXPECT_TRUE(testnewickequality(net, "(d,((f,(e)#i4)i5,(a,b)i1)i2,(c,#i4)i6)i3;"));
}

namespace {
void reroot_i1i2_and_compare(Network& net) {
    int idx = -1;
    for (size_t k = 0; k < net.edge.size(); ++k) {
        const Edge& e = net.E(net.edge[k]);
        std::vector<std::string> nm{net.N(e.nd[0]).name, net.N(e.nd[1]).name};
        std::sort(nm.begin(), nm.end());
        if (nm == std::vector<std::string>{"i1", "i2"}) {
            idx = static_cast<int>(k);
            break;
        }
    }
    ASSERT_GE(idx, 0);
    rootonedge_index(net, idx);
    semidirectnetwork(net);
    Network ref = reload_labelled_net(true);
    EXPECT_EQ(writenewick(net), writenewick(ref));
}
}  // namespace

TEST(rSPR, Reversible) {
    {
        Network net = reload_labelled_net(true);
        auto a = get_nodes(net, {"i1", "i3", "i6", "i2", "i4", "i5"});
        performrSPR(net, a[0], a[1], a[2], a[3], a[4], a[5]);
        auto b = get_nodes(net, {"i1", "i3", "i6", "i2", "i4", "i5"});
        performrSPR(net, b[0], b[4], b[5], b[3], b[1], b[2]);
        reroot_i1i2_and_compare(net);
    }
    {
        Network net = reload_labelled_net(true);
        auto a = get_nodes(net, {"i6", "i3", "i5", "i4", "i2", "i1"});
        performrSPR(net, a[0], a[1], a[2], a[3], a[4], a[5]);
        auto b = get_nodes(net, {"i6", "i3", "i5", "i4", "i2", "i1"});
        performrSPR(net, b[0], b[4], b[5], b[3], b[1], b[2]);
        reroot_i1i2_and_compare(net);
    }
}

TEST(rSPR, ManyRandomMoves) {
    Rng rng(0);
    Network net = reload_labelled_net(true);
    for (int j = 0; j < 1000; ++j) {
        std::string prev = writenewick(net);
        MoveParams p = samplerSPRparameters(net, rng);
        ASSERT_TRUE(p.valid);
        performrSPR(net, p.p[0].id, p.p[1].id, p.p[2].id, p.p[3].id, p.p[4].id, p.p[5].id);
        std::string now = writenewick(net);
        for (const auto& tip : net.tiplabels()) EXPECT_NE(now.find(tip), std::string::npos);
        EXPECT_NE(prev, now);
    }
}

// ------------------------------------------------------------------- rNNI

TEST(rNNI, VerifiedByHand) {
    Rng rng(0);
    {
        Network net = reload_labelled_net();
        auto n = get_nodes(net, {"i6", "i4", "i2", "i3"});
        performrNNI(net, n[0], n[1], n[2], n[3], 1, rng);
        EXPECT_TRUE(testnewickequality(net, "(((e,f)i5)#i4,(d,(c,#i4)i6)i3,(a,b)i1)i2;"));
    }
    {
        Network net = reload_labelled_net();
        auto n = get_nodes(net, {"i2", "i3", "i6", "i4"});
        performrNNI(net, n[0], n[1], n[2], n[3], 2, rng);
        EXPECT_TRUE(testnewickequality(net, "(((e,f)i5)#i4,(d,(c,#i4)i6)i3,(a,b)i1)i2;"));
    }
    {
        Network net = reload_labelled_net();
        auto n = get_nodes(net, {"i6", "f", "i4", "i5"});
        performrNNI(net, n[0], n[1], n[2], n[3], 3, rng);
        EXPECT_TRUE(testnewickequality(net, "((c,#H5)i6,(d,((e)#H5,f)i4)i3,(a,b)i1)i2;"));
    }
    {
        Network net = reload_labelled_net();
        auto n = get_nodes(net, {"d", "i6", "i3", "i4"});
        performrNNI(net, n[0], n[1], n[2], n[3], 4, rng);
        EXPECT_TRUE(testnewickequality(net, "((c,#H3)i6,(((e,f)i5,d)i4)#H3,(a,b)i1)i2;"));
    }
    {
        Network net = reload_labelled_net();
        auto n = get_nodes(net, {"c", "i3", "i6", "i4"});
        performrNNI(net, n[0], n[1], n[2], n[3], 4, rng);
        auto m = get_nodes(net, {"i1", "i3", "i2", "H6"});
        EXPECT_TRUE(isvalidrNNI4(net, m[0], m[1], m[2], m[3]));
    }
}

TEST(rNNI, Reversible) {
    Rng rng(0);
    auto twice = [&](std::initializer_list<const char*> fwd, std::initializer_list<const char*> back, int type) {
        Network net0 = reload_labelled_net(), net1 = reload_labelled_net(), net2 = reload_labelled_net();
        for (int j = 0; j < 100; ++j) {
            auto a = get_nodes(net1, fwd);
            performrNNI(net1, a[0], a[1], a[2], a[3], type, rng);
            EXPECT_NE(writenewick(net1), writenewick(net2));
            auto b = get_nodes(net2, fwd);
            performrNNI(net2, b[0], b[1], b[2], b[3], type, rng);
            EXPECT_EQ(writenewick(net1), writenewick(net2));
            auto c = get_nodes(net1, back);
            performrNNI(net1, c[0], c[1], c[2], c[3], type, rng);
            EXPECT_NE(writenewick(net1), writenewick(net2));
            auto d = get_nodes(net2, back);
            performrNNI(net2, d[0], d[1], d[2], d[3], type, rng);
            EXPECT_EQ(writenewick(net1), writenewick(net2));
            EXPECT_EQ(writenewick(net1), writenewick(net0));
        }
    };
    twice({"d", "i5", "i3", "i4"}, {"i5", "d", "i3", "i4"}, 1);
    twice({"i2", "i3", "i6", "i4"}, {"i3", "i2", "i6", "i4"}, 2);

    // (3) then (4)
    {
        Network net0 = reload_labelled_net(), net1 = reload_labelled_net(), net2 = reload_labelled_net();
        for (int j = 0; j < 100; ++j) {
            auto a = get_nodes(net1, {"i6", "f", "i4", "i5"});
            performrNNI(net1, a[0], a[1], a[2], a[3], 3, rng);
            EXPECT_NE(writenewick(net1), writenewick(net2));
            auto b = get_nodes(net2, {"i6", "f", "i4", "i5"});
            performrNNI(net2, b[0], b[1], b[2], b[3], 3, rng);
            EXPECT_EQ(writenewick(net1), writenewick(net2));
            net1.N(a[0]).name = "i6";
            net2.N(b[0]).name = "i6";
            net1.N(a[3]).name = "i5";
            net2.N(b[3]).name = "i5";
            performrNNI(net1, a[1], a[0], a[2], a[3], 4, rng);
            EXPECT_NE(writenewick(net1), writenewick(net2));
            performrNNI(net2, b[1], b[0], b[2], b[3], 4, rng);
            net1.N(a[2]).name = "i4";
            net2.N(b[2]).name = "i4";
            EXPECT_EQ(writenewick(net1), writenewick(net2));
            EXPECT_EQ(writenewick(net1), writenewick(net0));
        }
    }
    // (4) then (3)
    {
        Network net0 = reload_labelled_net(), net1 = reload_labelled_net(), net2 = reload_labelled_net();
        for (int j = 0; j < 100; ++j) {
            auto a = get_nodes(net1, {"d", "i6", "i3", "i4"});
            performrNNI(net1, a[0], a[1], a[2], a[3], 4, rng);
            EXPECT_NE(writenewick(net1), writenewick(net2));
            auto b = get_nodes(net2, {"d", "i6", "i3", "i4"});
            performrNNI(net2, b[0], b[1], b[2], b[3], 4, rng);
            EXPECT_EQ(writenewick(net1), writenewick(net2));
            net1.N(a[2]).name = "i3";
            net2.N(b[2]).name = "i3";
            auto c = get_nodes(net1, {"i6", "d", "i3", "i4"});
            performrNNI(net1, c[0], c[1], c[2], c[3], 3, rng);
            EXPECT_NE(writenewick(net1), writenewick(net2));
            auto d = get_nodes(net2, {"i6", "d", "i3", "i4"});
            performrNNI(net2, d[0], d[1], d[2], d[3], 3, rng);
            net1.N(a[3]).name = "i4";
            net2.N(b[3]).name = "i4";
            EXPECT_EQ(writenewick(net1), writenewick(net2));
            EXPECT_EQ(writenewick(net1), writenewick(net0));
        }
    }
}

TEST(rNNI, ManyMovesOfEachType) {
    Rng rng(1);
    Network net = reload_labelled_net();
    for (int type = 1; type <= 4; ++type) {
        for (int j = 0; j < 100; ++j) {
            std::string prev = writenewick(net);
            auto valid = allvalidrNNInodes(net, type);
            if (valid.empty()) {
                net = reload_labelled_net();
                valid = allvalidrNNInodes(net, type);
                prev = "";
            }
            ASSERT_FALSE(valid.empty());
            auto s = valid[rng.index(valid.size())];
            EXPECT_NO_THROW(performrNNI(net, s[0], s[1], s[2], s[3], type, rng));
            EXPECT_NO_THROW(writenewick(net));
            EXPECT_NE(prev, writenewick(net));
            bool allleaf = true;
            for (NodeId h : net.hybrid) allleaf &= net.N(net.getchild_node(h)).leaf;
            if (type == 3 && allleaf) net = reload_labelled_net();
        }
    }
}

TEST(rNNI, RandomMoves) {
    // performrandomrNNI!: pick a type with equal probabilities, falling back
    // to the other types when the chosen one has no valid move
    Rng rng(0);
    Network net = reload_labelled_net();
    for (int j = 1; j <= 10000; ++j) {
        if (j % 100 == 0) net = reload_labelled_net();
        std::string prev = writenewick(net);
        std::vector<int> types{1, 2, 3, 4};
        bool done = false;
        while (!done && !types.empty()) {
            size_t k = rng.index(types.size());
            int type = types[k];
            auto valid = allvalidrNNInodes(net, type);
            if (valid.empty()) {
                types.erase(types.begin() + static_cast<long>(k));
                continue;
            }
            auto s = valid[rng.index(valid.size())];
            performrNNI(net, s[0], s[1], s[2], s[3], type, rng);
            done = true;
        }
        ASSERT_TRUE(done);
        EXPECT_NO_THROW(writenewick(net));
        EXPECT_NE(prev, writenewick(net));
    }
}

TEST(rNNI, PastBugs) {
    {
        Network net = nw("(b:1.098,((((f:0.551,(e:0.276)#H1:0.276::0.75)i5:1.146,(c1:2.109,c2:1.557)i6:0.037)i7:0.926,((d1:1.602,d2:1.602)i10:0.086,a:1.687)i12:1.385)i13:1.585,#H1:0.0::0.25)i5554:1.585)i15;");
        semidirectnetwork(net);
        auto n = get_nodes(net, {"i13", "i5", "i5554", "H1"});
        EXPECT_FALSE(isvalidrNNI4(net, n[0], n[1], n[2], n[3]));
    }
    {
        Network net = nw("(((d1:1.602,d2:1.602)i10:0.086,c2:1.687)i12:1.385,b:4.267,(((c1:2.109,a:1.557)i6:0.037,((e:0.551,f:0.551)i5:0.573)#H1:0.573::0.75)i7:0.463,#H1:0.0::0.25)i-6079:0.463)i13;");
        net.isrooted = false;
        EXPECT_NO_THROW(allvalidrNNInodes(net, 1));
    }
    {
        Network net = nw("((a:1.146,((c1:2.109,(e:0.551,(f:0.276)#H1:0.5::0.765)i5:0.924)i6:0.005,#H1:0.5::0.235)i8004:0.002)i7:0.0,((d1:1.602,c2:1.602)i10:0.947,d2:1.687)i12:0.0,b:4.267)i13;");
        net.isrooted = false;
        auto n = get_nodes(net, {"i6", "f", "i8004", "H1"});
        EXPECT_TRUE(isvalidrNNI1(net, n[0], n[1], n[2], n[3]));
    }
    {
        Network net = nw("(a:1.146,(c1:2.109,((e:0.551,(f:0.276)#H1:0.5::0.501)i5:0.425,#H1:0.5::0.499)i3841:0.425)i6:0.007,(c2:0.947,(b:0.947,(d1:1.602,d2:1.602)i10:0.947)i12:0.0)i13:0.0)i7;");
        net.isrooted = false;
        auto n = get_nodes(net, {"i6", "i5", "i3841", "H1"});
        EXPECT_FALSE(isvalidrNNI2(net, n[0], n[1], n[2], n[3]));
    }
    {
        Network net = nw("((d1:1.602,d2:1.602)i10:0.947,((a:1.146,((c1:2.109,(f:0.551,(e:0.276)#H1:0.5::0.682)i5:1.12)i6:0.0,#H1:0.5::0.318)i-8938:0.007)i7:0.0,c2:1.687)i12:0.0,b:4.267)i13;");
        net.isrooted = false;
        auto n = get_nodes(net, {"i6", "i5", "i-8938", "H1"});
        EXPECT_FALSE(isvalidrNNI4(net, n[0], n[1], n[2], n[3]));
    }
    {
        Network net = nw("(((b:0.5,(((#H1:0.227::0.491,(c2:0.5,c1:0.875)i3:0.47)i12:0.396,(((d2:0.5,d1:0.75)i7:0.471,((e:0.5,f:1.0)Hi4:0.227)#H1:0.227::0.509)i6216:0.239)#H3:0.239::0.681)i4:0.521,#H3:0.239::0.319)i10076:0.615)i15:0.5,#H2:0.5::0.5)i7554:0.5,(a:0.261)#H2:0.5::0.5)i16;");
        net.isrooted = false;
        EXPECT_NO_THROW(allvalidrNNInodes(net, 2));
    }
    {
        Network net = nw("((((c1:0.219,c2:1.277)i5355:0.773,((e:0.551,f:0.551)i5:0.48,(#H1:0.5::0.091,(d1:0.65,d2:0.5)i10443:0.052)i8641:0.846)i6:0.126)i7793:1.435,(b:0.441)#H1:0.5::0.909)i7:0.5,a:2.634)i13;");
        net.isrooted = false;
        auto n = get_nodes(net, {"i13", "i8641", "i7", "H1"});
        EXPECT_FALSE(isvalidrNNI2(net, n[0], n[1], n[2], n[3]));
    }
    {
        Network net = nw("(a:1.146,(b:1.687,((((e:0.776,f:0.855)i5:0.49,(c1:0.289,(c2:0.482)#H1:0.5::0.905)i13:0.978)i6:0.125,(d1:1.602,d2:1.301)i10:0.818)i1346:0.789,#H1:0.5::0.095)i10616:0.417)i12:0.5)i7;");
        net.isrooted = false;
        auto n = get_nodes(net, {"i12", "i13", "i10616", "H1"});
        EXPECT_FALSE(isvalidrNNI2(net, n[0], n[1], n[2], n[3]));
    }
}

// --------------------------------------------------------------- all moves

TEST(AllMoves, ThousandsOfRandomMoves) {
    for (int z = 1; z <= 500; ++z) {
        Rng rng(static_cast<uint64_t>(z));
        int ntaxa = static_cast<int>(rng.randint(5, 40));
        int nhyb = static_cast<int>(rng.randint(0, ntaxa));
        Network t = random_tree(ntaxa, rng);
        semidirectnetwork(t);
        std::string prev = writenewick(t, true);
        for (int miter = 0; miter < 250; ++miter) {
            Network t2 = t.deepcopynetwork();
            std::vector<AttemptedMove> none;
            CFMatrix emptyQ;
            std::vector<double> deltas;
            std::pair<Move, MoveParams> mp;
            ASSERT_NO_THROW(mp = generatemoveproposal(t2, {}, none, nhyb, 0.0, emptyQ, deltas, rng, 0.0))
                << "z=" << z << " iter=" << miter << " newick=" << prev;
            if (mp.first == Move::none) break;
            try {
                applymove(t2, mp.first, mp.second, rng);
                removedegree2nodes(t2);
                while (shrink2cycles(t2)) {
                }
                semidirectnetwork(t2);
                Network tc = t2.deepcopynetwork();
                std::string now = writenewick(tc, true);
                t = std::move(t2);
                EXPECT_TRUE(mp.first == Move::rNNI2 || now != prev) << move_name(mp.first);
                for (NodeId n : t.node) EXPECT_NE(now.find(t.N(n).name), std::string::npos);
                prev = now;
            } catch (const std::exception& e) {
                FAIL() << "z=" << z << " iter=" << miter << " move=" << move_name(mp.first) << " error: " << e.what()
                       << "\nnewick before move: " << prev;
            }
        }
    }
}
