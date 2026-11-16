// Port of SNaQ.jl test/test_gradient_opt/test_opt_API.jl
#include <cmath>
#include <sstream>

#include "optimize.hpp"
#include "test_support.hpp"

using namespace snaqpp;
using namespace snaqpp::test;

namespace {
std::pair<Network, CFMatrix> get_data(double L) {
    Network net = read_trees("n1.netfile").at(0);
    semidirectnetwork(net);
    for (EdgeId e : net.edge) net.E(e).length = net.E(e).length > 0 ? L : 0;
    CFMatrix q = computeexpectedCFmatrix(net);
    return {std::move(net), std::move(q)};
}
}  // namespace

TEST(OptAPI, FitGetsCloseToTruth) {
    for (double L : {0.1, 0.25, 0.5, 1.0, 2.0, 5.0}) {
        for (int seed = 1; seed <= 10; ++seed) {
            auto [net, q] = get_data(L);
            OptInfo info = gatheroptimizationinfo(net, true);
            std::vector<double> params = gatherparams(net, info.param_map);
            Rng rng(static_cast<uint64_t>(seed));
            Network opt_net = net.clone();
            for (EdgeId e : opt_net.edge)
                if (opt_net.E(e).length != 0.0) opt_net.E(e).length = rng.rand() * L;
            for (NodeId h : opt_net.hybrid) {
                double g = rng.rand();
                g = std::max(g, 1.0 - g);
                opt_net.E(opt_net.getparentedge(h)).gamma = g;
                opt_net.E(opt_net.getparentedgeminor(h)).gamma = 1.0 - g;
            }
            FitOptions fo;
            fo.maxeval = 100000;
            Network copy = opt_net.deepcopynetwork();
            double first = fitnumericalparameters_all(copy, q, 0.0, fo);
            EquationSet es = findquartetequations(opt_net);
            double second = fitnumericalparameters(opt_net, es.eqns, q, 0.0, fo);
            OptInfo oi = gatheroptimizationinfo(opt_net, true);
            std::vector<double> opt_params = gatherparams(opt_net, oi.param_map);
            EXPECT_NEAR(first, second, 1e-5);
            double meanabs = 0.0;
            for (size_t k = 0; k < params.size(); ++k) meanabs += std::abs(opt_params[k] - params[k]);
            meanabs /= static_cast<double>(params.size());
            EXPECT_LT(meanabs, 1.0) << "L=" << L << " seed=" << seed;
        }
    }
}

namespace {
void strictly_improves(bool use_rho) {
    auto lines = read_lines("opt_improve.txt");
    ASSERT_EQ(lines.size(), 100u);
    for (size_t k = 0; k < lines.size(); ++k) {
        std::istringstream ss(lines[k]);
        double rho;
        std::string newick;
        ss >> rho >> newick;
        if (!use_rho) rho = 0.0;
        Network net = nw(newick);
        CFMatrix q = read_cf("opt_improve_" + std::to_string(k + 1) + ".cf");
        double before = computeSNaQscore_net(net, q, rho);
        EquationSet es = findquartetequations(net);
        FitOptions fo;
        fo.maxeval = 10;
        fitnumericalparameters(net, es.eqns, q, rho, fo);
        double after = computeSNaQscore_net(net, q, rho);
        EXPECT_GT(after, before) << "case " << k + 1;
    }
}
}  // namespace

TEST(OptAPI, LogPLStrictlyImproves) { strictly_improves(false); }
TEST(OptAPI, LogPLStrictlyImprovesWithRho) { strictly_improves(true); }
