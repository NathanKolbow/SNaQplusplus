// Port of SNaQ.jl test/test_gradient_opt/test_CF_recursive_blocks.jl
#include <cmath>

#include "test_support.hpp"

using namespace snaqpp;
using namespace snaqpp::test;

namespace {
RecursiveCFEquation make_eqn(bool can_coal, std::initializer_list<int> coal, int which, int divH) {
    RecursiveCFEquation e;
    e.can_coalesce_here = can_coal;
    for (int c : coal) e.coal_edges.push_back(c);
    e.which_coal = which;
    e.division_H = divH;
    return e;
}
}  // namespace

TEST(RecursiveCFEquation, Construction) {
    RecursiveCFEquation eqn = make_eqn(true, {0, 1, 2}, 1, -1);
    EXPECT_TRUE(eqn.can_coalesce_here);
    ASSERT_EQ(eqn.coal_edges.size(), 3u);
    EXPECT_EQ(eqn.coal_edges[0], 0);
    EXPECT_EQ(eqn.coal_edges[2], 2);
    EXPECT_EQ(eqn.which_coal, 1);
    EXPECT_EQ(eqn.division_H, -1);
    EXPECT_TRUE(eqn.divisions.empty());
}

TEST(RecursiveCFEquation, QuartetDataConstruction) {
    QuartetDataImpl q;
    q.eqn = make_eqn(true, {0, 1}, 1, -1);
    for (int p : {0, 1, 2}) q.relevant_params.push_back(p);
    q.q_taxa = {0, 1, 2, 3};
    EXPECT_EQ(q.relevant_params.size(), 3u);
    EXPECT_TRUE(q.relevant_params.contains(0));
    EXPECT_FALSE(q.relevant_params.contains(3));
}

TEST(RecursiveCFEquation, ExpectedCFTreelike) {
    QuartetDataImpl q;
    q.eqn = make_eqn(false, {0}, 1, -1);
    std::vector<double> params{0.5};
    auto e = computeexpectedCF(q, params, 0.0);
    EXPECT_GT(e[0], e[1]);
    EXPECT_NEAR(e[1], 1.0 / 3.0 * std::exp(-0.5), 1e-12);
    EXPECT_NEAR(e[0], 1 - 2.0 / 3.0 * std::exp(-0.5), 1e-12);
}

TEST(RecursiveCFEquation, LossAndGradient) {
    Network net = nw("((A,(B)#H1),((C,D), #H1));");
    for (EdgeId e : net.edge) {
        net.E(e).length = 1.0;
        net.E(e).gamma = !net.E(e).hybrid ? -1 : 0.5;
    }
    semidirectnetwork(net);
    CFMatrix q;
    q.v = {0.6, 0.2, 0.2};
    EquationSet es = findquartetequations(net);
    std::vector<double> grad;
    double loss = computelossandgradient(es.eqns, es.info.params, grad, q, rhotoalpha(0.0));
    EXPECT_NEAR(loss, -0.20561335094444322, 1e-12);
    const double eps = 1e-6;
    for (size_t i = 0; i < es.info.params.size(); ++i) {
        auto pp = es.info.params, pm = es.info.params;
        pp[i] += eps;
        pm[i] -= eps;
        double num = (computeSNaQscore(es.eqns, pp, q, 0.0) - computeSNaQscore(es.eqns, pm, q, 0.0)) / (2 * eps);
        if (std::abs(num) > 1e-6) EXPECT_NEAR(grad[i], num, 0.1 * std::abs(num));
    }
}

TEST(RecursiveCFEquation, RecursionSimpleCase) {
    QuartetDataImpl q;
    q.eqn = make_eqn(false, {0}, 1, -1);
    // gradient of eCF1 and eCF2 with respect to the single branch, checked
    // through finite differences of computeexpectedCF
    std::vector<double> p{0.5};
    auto e = computeexpectedCF(q, p, 0.0);
    EXPECT_NEAR(e[0], 1 - 2.0 / 3.0 * std::exp(-0.5), 1e-12);
    EXPECT_NEAR(e[1], 1.0 / 3.0 * std::exp(-0.5), 1e-12);
    const double h = 1e-7;
    std::vector<double> pp{0.5 + h}, pm{0.5 - h};
    auto ep = computeexpectedCF(q, pp, 0.0), em = computeexpectedCF(q, pm, 0.0);
    EXPECT_NEAR((ep[0] - em[0]) / (2 * h), 2.0 / 3.0 * std::exp(-0.5), 1e-6);
    EXPECT_NEAR((ep[1] - em[1]) / (2 * h), -1.0 / 3.0 * std::exp(-0.5), 1e-6);
    double e3p = 1 - ep[0] - ep[1], e3m = 1 - em[0] - em[1];
    EXPECT_NEAR((e3p - e3m) / (2 * h), -1.0 / 3.0 * std::exp(-0.5), 1e-6);
}
