#include "optimize.hpp"

#include "classes.hpp"

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

namespace {

struct ClassObjective {
    const ClassBatch* batch;
    const std::vector<int>* pidof;
    std::vector<double> params, grad;  // indexed by pid
    double alpha;
    LossWorkspace ws;
    int nthreads;
    int evals = 0;
};

double classobjective(unsigned n, const double* x, double* grad, void* data) {
    auto* d = static_cast<ClassObjective*>(data);
    ++d->evals;
    const std::vector<int>& pidof = *d->pidof;
    for (unsigned j = 0; j < n; ++j) d->params[static_cast<size_t>(pidof[j])] = x[j];
    const double f = evalclasses(*d->batch, d->params.data(), d->grad.data(), d->alpha, d->ws, d->nthreads);
    if (grad)
        for (unsigned j = 0; j < n; ++j) grad[j] = d->grad[static_cast<size_t>(pidof[j])];
    return f;
}

// Solves (B + lambda (diag(B) + floor)) d = g on the free variables by
// Cholesky factorization. `B` (n x n) is symmetric positive semidefinite.
// Returns false if the matrix is not positive definite.
bool dampedsolve(const std::vector<double>& B, const std::vector<double>& g, const std::vector<int>& freeidx,
                 double lambda, double floor, std::vector<double>& L, std::vector<double>& d, size_t n) {
    const size_t m = freeidx.size();
    L.assign(m * m, 0.0);
    for (size_t a = 0; a < m; ++a) {
        const size_t ia = static_cast<size_t>(freeidx[a]);
        for (size_t b = 0; b <= a; ++b) L[a * m + b] = B[ia * n + static_cast<size_t>(freeidx[b])];
        L[a * m + a] += lambda * (B[ia * n + ia] + floor);
    }
    for (size_t a = 0; a < m; ++a) {
        for (size_t b = 0; b <= a; ++b) {
            double sum = L[a * m + b];
            for (size_t k = 0; k < b; ++k) sum -= L[a * m + k] * L[b * m + k];
            if (a == b) {
                if (!(sum > 0.0)) return false;
                L[a * m + a] = std::sqrt(sum);
            } else {
                L[a * m + b] = sum / L[b * m + b];
            }
        }
    }
    d.assign(m, 0.0);
    for (size_t a = 0; a < m; ++a) {  // L y = g
        double sum = g[static_cast<size_t>(freeidx[a])];
        for (size_t k = 0; k < a; ++k) sum -= L[a * m + k] * d[k];
        d[a] = sum / L[a * m + a];
    }
    for (size_t a = m; a-- > 0;) {  // L' d = y
        double sum = d[a];
        for (size_t k = a + 1; k < m; ++k) sum -= L[k * m + a] * d[k];
        d[a] = sum / L[a * m + a];
    }
    return true;
}

// Projected Newton method with Levenberg-Marquardt damping, for the
// maximization of the class likelihood within the parameters' bounds.
//
// The curvature matrix B of evalclasses is exact for tree-like classes, whose
// terms are concave in the branch lengths, and a Gauss-Newton approximation
// for reticulate classes; it is positive semidefinite, so every damped step
// is an ascent direction. Variables at a bound whose gradient points outside
// are held fixed; a step is accepted only if it increases the score.
double newtonfit(const ClassBatch& batch, const std::vector<int>& pidof, std::vector<double>& params, double alpha,
                 std::vector<double>& x, const OptInfo& info, const FitOptions& fo, int maxeval, int& evals,
                 bool& converged) {
    const size_t n = x.size();
    std::vector<int> dense(static_cast<size_t>(batch.npid), -1);
    for (size_t j = 0; j < n; ++j) dense[static_cast<size_t>(pidof[j])] = static_cast<int>(j);
    std::vector<double> gpid(static_cast<size_t>(batch.npid)), g(n), B(n * n), xt(n), L, d;
    LossWorkspace ws;
    // The score at a point; then, only if the point is kept, its gradient and
    // curvature (a rejected trial step needs neither).
    auto score = [&](const std::vector<double>& at) {
        for (size_t j = 0; j < n; ++j) params[static_cast<size_t>(pidof[j])] = at[j];
        ++evals;
        return evalclassvalues(batch, params.data(), alpha, ws, 1);
    };
    auto derivatives = [&](std::vector<double>& grad, std::vector<double>& curv) {  // at the last point scored
        evalclassderivatives(batch, params.data(), gpid.data(), alpha, ws, 1, dense.data(), static_cast<int>(n),
                             curv.data());
        for (size_t j = 0; j < n; ++j) grad[j] = gpid[static_cast<size_t>(pidof[j])];
    };
    double f = score(x);
    converged = false;
    if (n == 0 || !std::isfinite(f)) {
        converged = n == 0;
        return f;
    }
    double lambda = 1e-3;
    int smallgains = 0;
    std::vector<int> freeidx;
    std::vector<double> step(n);
    while (evals < maxeval) {
        derivatives(g, B);  // x is the last point scored
        // free variables: not held at a bound by the gradient
        freeidx.clear();
        double maxdiag = 0.0, pg = 0.0;
        for (size_t j = 0; j < n; ++j) {
            const bool atlow = x[j] <= info.LB[j] + 1e-10 && g[j] < 0.0;
            const bool athigh = x[j] >= info.UB[j] - 1e-10 && g[j] > 0.0;
            if (atlow || athigh) continue;
            freeidx.push_back(static_cast<int>(j));
            maxdiag = std::max(maxdiag, B[j * n + j]);
            pg = std::max(pg, std::abs(g[j]));
        }
        if (freeidx.empty() || pg <= 1e-11 * (1.0 + std::abs(f))) {
            converged = true;
            break;
        }
        const double floor = 1e-8 * maxdiag + 1e-12;
        bool accepted = false, stalled = false;
        double fnew = f, usedlambda = lambda;
        for (int attempt = 0; attempt < 30 && evals < maxeval; ++attempt) {
            if (!dampedsolve(B, g, freeidx, lambda, floor, L, d, n)) {
                lambda *= 10.0;
                continue;
            }
            // step, limited per parameter and kept within the bounds
            xt = x;
            std::fill(step.begin(), step.end(), 0.0);
            double maxstep = 0.0, linear = 0.0;
            for (size_t a = 0; a < freeidx.size(); ++a) {
                const size_t j = static_cast<size_t>(freeidx[a]);
                const double limit = info.idx_obj[j].is_hybrid_node ? 0.25 : 3.0;
                double v = x[j] + std::max(-limit, std::min(limit, d[a]));
                v = std::max(info.LB[j], std::min(info.UB[j], v));
                xt[j] = v;
                step[j] = v - x[j];
                maxstep = std::max(maxstep, std::abs(step[j]));
                linear += g[j] * step[j];
            }
            if (maxstep <= 1e-14 || linear <= 0.0) {
                stalled = true;  // no room to move
                break;
            }
            // gain predicted by the quadratic model
            double quad = 0.0;
            for (size_t a = 0; a < freeidx.size(); ++a) {
                const size_t j = static_cast<size_t>(freeidx[a]);
                if (step[j] == 0.0) continue;
                double row = 0.0;
                for (size_t c = 0; c < freeidx.size(); ++c) {
                    const size_t k = static_cast<size_t>(freeidx[c]);
                    row += B[j * n + k] * step[k];
                }
                quad += step[j] * row;
            }
            const double predicted = linear - 0.5 * quad;
            fnew = score(xt);
            usedlambda = lambda;
            const double ratio = (std::isfinite(fnew) && predicted > 0.0) ? (fnew - f) / predicted : -1.0;
            if (std::isfinite(fnew) && fnew > f) {
                accepted = true;
                if (ratio > 0.75)
                    lambda = std::max(lambda * 0.2, 1e-10);
                else if (ratio < 0.25)
                    lambda *= 4.0;
                break;
            }
            lambda *= 8.0;
            if (lambda > 1e14) {
                stalled = true;
                break;
            }
        }
        if (!accepted) {
            // no damped step improves the score: x is a (local) maximum to working precision,
            // unless the evaluations ran out
            converged = stalled || evals < maxeval;
            break;
        }
        const double gain = fnew - f;
        x.swap(xt);
        f = fnew;
        // A small gain means convergence only if the step was close to a full
        // Newton step; heavily damped steps are small whatever the distance to
        // the optimum (but several of them in a row are not worth continuing).
        if (gain <= fo.ftolAbs || gain <= fo.ftolRel * std::abs(f)) {
            ++smallgains;
            if (usedlambda <= 1e-2 || smallgains >= 8) {
                converged = usedlambda <= 1e-2;
                break;
            }
        } else {
            smallgains = 0;
        }
    }
    return f;
}

}  // namespace

double fitclasses(Network& net, const QuartetClasses& cls, double rho, const FitOptions& fo, int nthreads,
                  FitInfo* fitinfo) {

    if (!(rho >= 0 && rho <= 1)) throw SnaqError("ρ must be in range [0, 1]");
    if (cls.npid != net.next_pid) throw SnaqError("fitclasses: the classes were built for another network");
    if (!cls.complete()) throw SnaqError("fitclasses: finishclasses must be called first");
    for (EdgeId e : net.edge) {
        Edge& ed = net.E(e);
        if (std::isnan(ed.length)) ed.length = 0.0;
        if (net.N(net.getchild(e)).leaf) continue;
        ed.length = std::max(ed.length, 1e-5);  // starting optimization on a boundary can lead to failure
    }
    OptInfo info = gatheroptimizationinfo(net, false);
    const unsigned narg = static_cast<unsigned>(info.narg);
    const bool newton = fo.method == FitMethod::newton;
    const int maxeval = fo.maxeval >= 0 ? fo.maxeval : (newton ? 200 : std::max(30, net.numtaxa));
    std::vector<int> pidof(narg);
    for (unsigned j = 0; j < narg; ++j) {
        const ParamObj& o = info.idx_obj[j];
        pidof[j] = o.is_hybrid_node ? net.N(o.id).pid : net.E(o.id).pid;
    }
    std::vector<double> x(narg);
    for (unsigned j = 0; j < narg; ++j) {
        double v = std::min(info.UB[j], info.params[j]);
        v = std::min(info.UB[j] - 1e-12, v);
        v = std::max(info.LB[j] + 1e-12, v);
        x[j] = v;
    }
    ClassObjective data{&cls.batch, &pidof, pidparams(net, cls.npid),
                        std::vector<double>(static_cast<size_t>(cls.npid)), rhotoalpha(rho), {}, nthreads};
    double maxf = std::numeric_limits<double>::quiet_NaN();
    bool converged = false;
    auto lbfgs = [&](int evalsleft) {
        NloptHandle h(narg);
        nlopt_set_maxeval(h.o, evalsleft);
        nlopt_set_ftol_rel(h.o, fo.ftolRel);
        nlopt_set_ftol_abs(h.o, fo.ftolAbs);
        nlopt_set_xtol_rel(h.o, fo.xtolRel);
        nlopt_set_xtol_abs1(h.o, fo.xtolAbs);
        if (narg > 0) {
            nlopt_set_initial_step(h.o, info.init_steps.data());
            nlopt_set_lower_bounds(h.o, info.LB.data());
            nlopt_set_upper_bounds(h.o, info.UB.data());
        }
        nlopt_set_max_objective(h.o, classobjective, &data);
        // NLopt.jl ignores the return code (except forced stops from callbacks).
        const nlopt_result ret = nlopt_optimize(h.o, x.data(), &maxf);
        converged = ret > 0 && ret != NLOPT_MAXEVAL_REACHED;
    };
    if (newton) {
        maxf = newtonfit(cls.batch, pidof, data.params, data.alpha, x, info, fo, maxeval, data.evals, converged);
        // The curvature of reticulate classes is approximate (Gauss-Newton), and
        // Newton steps can stall on it far from the optimum: L-BFGS finishes
        // from where they stopped, and stops at once if they had converged.
        // Trees need no such check: their curvature is exact and the score concave.
        if (!cls.batch.gen.empty() && narg > 0 && std::isfinite(maxf) && data.evals < maxeval) {
            const double fnewton = maxf;
            const std::vector<double> xnewton = x;
            const bool newtonconverged = converged;
            for (unsigned j = 0; j < narg; ++j) {  // NLopt needs a start strictly inside the bounds
                x[j] = std::min(info.UB[j] - 1e-12, std::max(info.LB[j] + 1e-12, x[j]));
            }
            lbfgs(maxeval - data.evals);
            if (!(maxf >= fnewton)) {  // (also if NLopt failed and left maxf undefined)
                maxf = fnewton;
                x = xnewton;
            }
            // L-BFGS started at an optimum can stop on an error (its line search
            // finds nothing better, to rounding): the fit has converged all the same.
            if (!converged && newtonconverged && data.evals < maxeval &&
                maxf - fnewton <= std::max(fo.ftolAbs, fo.ftolRel * std::abs(fnewton)))
                converged = true;
        }
    } else {
        lbfgs(maxeval);
    }
    setX(net, x.data(), info.idx_obj);
    if (maxf == -std::numeric_limits<double>::infinity()) throw SnaqError("Optimization error: maxf == -Inf");

    bool anyflip = false;
    for (NodeId hyb : net.hybrid) {
        auto par = net.getparents(hyb);
        EdgeId e1 = net.getconnectingedge(hyb, par[0]);
        EdgeId e2 = net.getconnectingedge(hyb, par[1]);
        const bool was = net.E(e1).ismajor;
        if (net.E(e1).gamma > 0.5) {
            net.E(e1).ismajor = true;
            net.E(e2).ismajor = false;
        } else if (net.E(e1).gamma < 0.5) {
            net.E(e1).ismajor = false;
            net.E(e2).ismajor = true;
        }
        anyflip |= was != net.E(e1).ismajor;
    }
    if (fitinfo) {
        fitinfo->flipped = anyflip;
        fitinfo->converged = converged;
        fitinfo->evals = data.evals;
    }
    net.fscore = maxf;
    return maxf;
}

double fitnumericalparameters_all(Network& net, const CFMatrix& obs, double rho, const FitOptions& opt,
                                  int nthreads) {
    fixnegativeedges(net);
    checkgammas(net);
    ClassOptions co;
    co.qlogq = obsqlogq(obs);
    co.nthreads = nthreads;
    QuartetClasses cls = buildclasses(net, QuartetList::allquartets(net.numtaxa), obs, co);
    return fitclasses(net, cls, rho, opt, nthreads);
}

}  // namespace csnaq
