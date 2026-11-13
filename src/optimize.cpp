#include "optimize.hpp"

#include <nlopt.h>

#include <algorithm>
#include <cmath>
#include <limits>

namespace csnaq {

namespace {

struct ObjectiveData {
    Network* net;
    const OptInfo* info;
    const TreeQuartetBatch* batch;
    const std::vector<QuartetData>* eqns;
    const CFMatrix* obs;
    double alpha;
    LossWorkspace* ws;
    int nthreads;
    std::vector<double> gradbuf;
};

double objective(unsigned n, const double* x, double* grad, void* data) {
    auto* d = static_cast<ObjectiveData*>(data);
    setX(*d->net, x, d->info->idx_obj);
    double* g = grad ? grad : d->gradbuf.data();
    return computelossandgradient(*d->batch, *d->eqns, x, n, g, *d->obs, d->alpha, *d->ws, d->nthreads);
}

struct NloptHandle {
    nlopt_opt o;
    explicit NloptHandle(unsigned n) : o(nlopt_create(NLOPT_LD_LBFGS, n)) {
        if (!o) throw SnaqError("could not create NLopt optimizer");
    }
    ~NloptHandle() { nlopt_destroy(o); }
};

}  // namespace

double fitnumericalparameters(Network& net, const std::vector<QuartetData>& eqns, const CFMatrix& obs, double rho,
                              const FitOptions& fo, int nthreads) {
    if (!(rho >= 0 && rho <= 1)) throw SnaqError("ρ must be in range [0, 1]");
    for (EdgeId e : net.edge) {
        Edge& ed = net.E(e);
        if (std::isnan(ed.length)) ed.length = 0.0;
        if (net.N(net.getchild(e)).leaf) continue;
        ed.length = std::max(ed.length, 1e-5);  // starting optimization on a boundary can lead to failure
    }
    const double alpha = rhotoalpha(rho);
    OptInfo info = gatheroptimizationinfo(net, false);
    const unsigned narg = static_cast<unsigned>(info.narg);
    const int maxeval = fo.maxeval >= 0 ? fo.maxeval : std::max(30, net.numtaxa);

    NloptHandle h(narg);
    nlopt_set_maxeval(h.o, maxeval);
    nlopt_set_ftol_rel(h.o, fo.ftolRel);
    nlopt_set_ftol_abs(h.o, fo.ftolAbs);
    nlopt_set_xtol_rel(h.o, fo.xtolRel);
    nlopt_set_xtol_abs1(h.o, fo.xtolAbs);
    if (narg > 0) {
        nlopt_set_initial_step(h.o, info.init_steps.data());
        nlopt_set_lower_bounds(h.o, info.LB.data());
        nlopt_set_upper_bounds(h.o, info.UB.data());
    }

    std::vector<double> x(narg);
    for (unsigned j = 0; j < narg; ++j) {
        double v = std::min(info.UB[j], info.params[j]);
        v = std::min(info.UB[j] - 1e-12, v);
        v = std::max(info.LB[j] + 1e-12, v);
        x[j] = v;
    }
    TreeQuartetBatch batch = treequartetbatch(eqns, obs);
    LossWorkspace ws;
    ObjectiveData data{&net, &info, &batch, &eqns, &obs, alpha, &ws, nthreads, std::vector<double>(narg + 1)};
    nlopt_set_max_objective(h.o, objective, &data);
    double maxf = std::numeric_limits<double>::quiet_NaN();
    // NLopt.jl ignores the return code (except forced stops from callbacks).
    nlopt_optimize(h.o, x.data(), &maxf);

    setX(net, x.data(), info.idx_obj);
    if (maxf == -std::numeric_limits<double>::infinity()) throw SnaqError("Optimization error: maxf == -Inf");

    for (NodeId hyb : net.hybrid) {
        auto par = net.getparents(hyb);
        EdgeId e1 = net.getconnectingedge(hyb, par[0]);
        EdgeId e2 = net.getconnectingedge(hyb, par[1]);
        if (net.E(e1).gamma > 0.5) {
            net.E(e1).ismajor = true;
            net.E(e2).ismajor = false;
        } else if (net.E(e1).gamma < 0.5) {
            net.E(e1).ismajor = false;
            net.E(e2).ismajor = true;
        }
    }
    net.fscore = maxf;
    return maxf;
}

double fitnumericalparameters_all(Network& net, const CFMatrix& obs, double rho, const FitOptions& opt,
                                  int nthreads) {
    EquationSet es = findquartetequations(net, nthreads);
    return fitnumericalparameters(net, es.eqns, obs, rho, opt, nthreads);
}

}  // namespace csnaq
