#include "likelihood.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <unordered_map>

#include "classes.hpp"
#include "parallel.hpp"
#include "snaq_network.hpp"

namespace snaqpp {

double rhotoalpha(double rho) { return rho == 0.0 ? std::numeric_limits<double>::infinity() : (1.0 - rho) / rho; }
double alphatorho(double alpha) { return std::isinf(alpha) ? 0.0 : 1.0 / (1.0 + alpha); }

namespace {

inline double coalbranchsum(const RecursiveCFEquation& eqn, const double* params) {
    double s = 0.0;
    for (int p : eqn.coal_edges) s += params[p];
    return s;
}

}  // namespace

double quadsplitprobability(int type, double g, double a) {
    if (std::isinf(a)) {
        if (type == 1) return g * g;
        if (type == 2) return (1 - g) * (1 - g);
        return g * (1 - g);
    } else if (a == 0.0) {
        if (type == 1) return g;
        if (type == 2) return 1 - g;
        return 0;
    }
    if (type == 1) return g * (1 / (a + 1) + a / (a + 1) * g);
    if (type == 2) return (1 - g) * (1 / (a + 1) + a / (a + 1) * (1 - g));
    if (type == 3 || type == 4) return g * (a / (a + 1)) * (1 - g);
    throw SnaqError("Found impossible type");
}

double quadsplitprobabilitygradient(int type, double g, double a) {
    if (std::isinf(a)) {
        if (type == 1) return 2 * g;
        if (type == 2) return -2 * (1 - g);
        return 1 - 2 * g;
    } else if (a == 0.0) {
        if (type == 1) return 1;
        if (type == 2) return -1;
        return 0.0;
    }
    if (type == 1) return 1 / (a + 1) + 2 * g * a / (a + 1);
    if (type == 2) return -1 / (a + 1) - 2 * (1 - g) * a / (a + 1);
    if (type == 3 || type == 4) return (a / (a + 1)) * (1 - 2 * g);
    throw SnaqError("Found impossible type");
}

namespace {

// Expected CFs (e1, e2, e3) of an equation. If `gs` is given (nparam x 3,
// row-major), w * d e_k / d p is added to gs[3*p + k] for every parameter p of
// the equation.
//
// One pass: `w`, the weight of this sub-equation in the whole equation (the
// product of the probabilities on the way to it), goes down, and values come
// back up. A parameter's derivative is the sum, over the places where it
// appears, of the weight of the place times the local derivative there.
// (SNaQ.jl's running_gradient.jl obtains the same derivatives with a list of
// the parameters seen on the way down, updated at every division.)
std::array<double, 3> recur(const RecursiveCFEquation& eqn, const double* params, double* gs, double w, double alpha) {
    const auto& coal = eqn.coal_edges;
    const int wc = eqn.which_coal;
    if (eqn.division_H == -1) {
        const double E = std::exp(-coalbranchsum(eqn, params));
        const double res = 1 - 2.0 / 3.0 * E, oth = 1.0 / 3.0 * E;
        if (gs) {
            const double wo = w * oth;
            for (int p : coal) {
                gs[3 * p] += wc == 1 ? 2.0 * wo : -wo;
                gs[3 * p + 1] += wc == 2 ? 2.0 * wo : -wo;
                gs[3 * p + 2] += wc == 3 ? 2.0 * wo : -wo;
            }
        }
        if (wc == 1) return {res, oth, oth};
        if (wc == 2) return {oth, res, oth};
        return {oth, oth, res};
    } else if (eqn.divisions.size() == 4) {
        double early = std::exp(-coalbranchsum(eqn, params));
        early = std::max(early, 1e-9);
        const int H = eqn.division_H;
        const double g = params[H];
        const double we = w * early;
        double below[3] = {0.0, 0.0, 0.0};  // early * sum of split_prob * division
        double dg[3] = {0.0, 0.0, 0.0};     // sum of split_grad * division
        for (int di = 1; di <= 4; ++di) {
            const double split_grad = quadsplitprobabilitygradient(di, g, alpha);
            const double split_prob = quadsplitprobability(di, g, alpha);
            const auto r = recur(eqn.divisions[static_cast<size_t>(di - 1)], params, gs, we * split_prob, alpha);
            const double fct = early * split_prob;
            for (int k = 0; k < 3; ++k) {
                below[k] += fct * r[static_cast<size_t>(k)];
                dg[k] += split_grad * r[static_cast<size_t>(k)];
            }
        }
        const bool here = eqn.can_coalesce_here && wc >= 1;
        if (gs) {
            for (int k = 0; k < 3; ++k) gs[3 * H + k] += we * dg[k];
            for (int p : coal) {
                for (int k = 0; k < 3; ++k) gs[3 * p + k] -= w * below[k];
                if (here) gs[3 * p + wc - 1] += we;
            }
        }
        std::array<double, 3> e = {below[0], below[1], below[2]};
        if (here) e[static_cast<size_t>(wc - 1)] += 1 - early;
        return e;
    } else {
        if (eqn.can_coalesce_here) throw SnaqError("Can coalesce w/ 1 taxa splitting at a single hybrid??");
        const int H = eqn.division_H;
        double g = params[H];
        g = std::max(1e-9, g);
        g = std::min(1.0 - 1e-9, g);
        const auto r = recur(eqn.divisions[0], params, gs, w * g, alpha);
        const auto s = recur(eqn.divisions[1], params, gs, w * (1 - g), alpha);
        std::array<double, 3> e;
        for (size_t k = 0; k < 3; ++k) e[k] = g * r[k] + (1 - g) * s[k];
        if (gs)
            for (size_t k = 0; k < 3; ++k) gs[3 * static_cast<size_t>(H) + k] += w * (r[k] - s[k]);
        return e;
    }
}

// Loss + gradient of one quartet no reticulation divides (non-batch path).
double treequartetlossandgradient(const RecursiveCFEquation& eqn, const double* params, const double* q,
                                  double* grad) {
    double branchsum = 0.0;
    for (int p : eqn.coal_edges) branchsum += params[p];
    const double E = std::exp(-branchsum);
    const int wc = eqn.which_coal;
    double e1 = wc == 1 ? 1 - 2.0 / 3.0 * E : 1.0 / 3.0 * E;
    double e2 = wc == 2 ? 1 - 2.0 / 3.0 * E : 1.0 / 3.0 * E;
    double e3 = 1.0 - e1 - e2;
    e1 = std::max(e1, 1e-9);
    e2 = std::max(e2, 1e-9);
    e3 = std::max(e3, 1e-9);
    const double q1 = q[0], q2 = q[1], q3 = q[2];
    double loss = ((q1 > 0) ? q1 * std::log(e1 / q1) : 0.0) + ((q2 > 0) ? q2 * std::log(e2 / q2) : 0.0) +
                  ((q3 > 0) ? q3 * std::log(e3 / q3) : 0.0);
    double gc = q1 * (wc == 1 ? 2.0 / 3.0 * E : -1.0 / 3.0 * E) / e1 +
                q2 * (wc == 2 ? 2.0 / 3.0 * E : -1.0 / 3.0 * E) / e2 +
                q3 * (wc == 3 ? 2.0 / 3.0 * E : -1.0 / 3.0 * E) / e3;
    for (int p : eqn.coal_edges) grad[p] += gc;
    return loss;
}

struct ThreadScratch {
    std::vector<double> quartet_grad;  // nparam x 3
    std::vector<double> total_grad;    // nparam
    explicit ThreadScratch(int np)
        : quartet_grad(3 * static_cast<size_t>(np), 0.0), total_grad(static_cast<size_t>(np), 0.0) {}
};

// General quartet: loss, gradient accumulated into ts.total_grad.
double quartetlossandgradient_general(const QuartetDataImpl& qd, const double* params, const double* q, double alpha,
                                      ThreadScratch& ts) {
    double* ig = ts.quartet_grad.data();
    for (int p : qd.relevant_params) {
        ig[3 * p] = 0.0;
        ig[3 * p + 1] = 0.0;
        ig[3 * p + 2] = 0.0;
    }
    auto r = recur(qd.eqn, params, ig, 1.0, alpha);
    double e1 = r[0], e2 = r[1];
    double e3 = 1.0 - e1 - e2;
    e1 = std::max(e1, 1e-9);
    e2 = std::max(e2, 1e-9);
    e3 = std::max(e3, 1e-9);
    const double i1 = 1.0 / e1, i2 = 1.0 / e2, i3 = 1.0 / e3;
    const double q1 = q[0], q2 = q[1], q3 = q[2];
    double loss = ((q1 > 0) ? q1 * std::log(e1 / q1) : 0.0) + ((q2 > 0) ? q2 * std::log(e2 / q2) : 0.0) +
                  ((q3 > 0) ? q3 * std::log(e3 / q3) : 0.0);
    for (int p : qd.relevant_params)
        ts.total_grad[static_cast<size_t>(p)] += q1 * ig[3 * p] * i1 + q2 * ig[3 * p + 1] * i2 + q3 * ig[3 * p + 2] * i3;
    return loss;
}

}  // namespace

struct LossWorkspace::Impl {
    int nparam = -1;
    std::vector<ThreadScratch> scratch;
    std::vector<double> tree_grad;  // nthreads x nparam
    std::vector<double> partial;    // per-thread partial losses
    // evalclassvalues -> evalclassderivatives
    std::vector<double> treeE;      // per tree-like class: exp(-path length)
    EqnDag::Work dag;
    std::vector<double> j1, j2, u1, u2;  // curvature of one class
    std::vector<int> jidx;
};

LossWorkspace::~LossWorkspace() { delete impl_; }

LossWorkspace::Impl* LossWorkspace::get(int nparam, int nthreads) {
    if (!impl_) impl_ = new Impl();
    if (impl_->nparam != nparam || static_cast<int>(impl_->scratch.size()) < nthreads) {
        impl_->nparam = nparam;
        impl_->scratch.clear();
        for (int t = 0; t < std::max(1, nthreads); ++t) impl_->scratch.emplace_back(nparam);
        impl_->tree_grad.assign(static_cast<size_t>(std::max(1, nthreads)) * static_cast<size_t>(nparam), 0.0);
        impl_->partial.assign(static_cast<size_t>(std::max(1, nthreads)), 0.0);
    }
    return impl_;
}

namespace {

struct VecHash {
    size_t operator()(const std::vector<int>& v) const {
        uint64_t h = 1469598103934665603ULL;
        for (int x : v) {
            h ^= static_cast<uint64_t>(static_cast<uint32_t>(x));
            h *= 1099511628211ULL;
        }
        return static_cast<size_t>(h ^ (h >> 29));
    }
};

void serialize_eqn(const RecursiveCFEquation& e, std::vector<int>& out) {
    out.push_back(e.can_coalesce_here ? 1 : 0);
    out.push_back(e.which_coal);
    out.push_back(e.division_H);
    std::vector<int> c(e.coal_edges.begin(), e.coal_edges.end());
    std::sort(c.begin(), c.end());
    out.push_back(static_cast<int>(c.size()));
    out.insert(out.end(), c.begin(), c.end());
    out.push_back(static_cast<int>(e.divisions.size()));
    for (const auto& d : e.divisions) serialize_eqn(d, out);
}

void eqn_params(const RecursiveCFEquation& e, std::vector<int>& out) {
    for (int p : e.coal_edges)
        if (std::find(out.begin(), out.end(), p) == out.end()) out.push_back(p);
    if (e.division_H >= 0 && std::find(out.begin(), out.end(), e.division_H) == out.end())
        out.push_back(e.division_H);
    for (const auto& d : e.divisions) eqn_params(d, out);
}

inline double qlogq3(double q1, double q2, double q3) {
    return ((q1 > 0) ? q1 * std::log(q1) : 0.0) + ((q2 > 0) ? q2 * std::log(q2) : 0.0) +
           ((q3 > 0) ? q3 * std::log(q3) : 0.0);
}

}  // namespace

TreeQuartetBatch treequartetbatch(const std::vector<QuartetData>& qdata, const CFMatrix& obs) {
    TreeQuartetBatch b;
    b.nquartets = qdata.size();
    std::unordered_map<std::vector<int>, size_t, VecHash> tree_groups, gen_groups;
    std::vector<std::vector<int>> tree_keys;
    std::vector<int> key;
    for (size_t j = 0; j < qdata.size(); ++j) {
        const auto& eqn = qdata[j]->eqn;
        const double q1 = obs(j, 0), q2 = obs(j, 1), q3 = obs(j, 2);
        if (eqn.division_H == -1) {
            key.assign(eqn.coal_edges.begin(), eqn.coal_edges.end());
            std::sort(key.begin(), key.end());
            auto it = tree_groups.find(key);
            size_t g;
            if (it == tree_groups.end()) {
                g = b.qres.size();
                tree_groups.emplace(key, g);
                tree_keys.push_back(key);
                b.qres.push_back(0.0);
                b.qoth.push_back(0.0);
                b.qlogq.push_back(0.0);
            } else {
                g = it->second;
            }
            const int wc = eqn.which_coal;
            const double r = wc == 1 ? q1 : (wc == 2 ? q2 : q3);
            b.qres[g] += r;
            b.qoth[g] += (q1 + q2 + q3) - r;
            b.qlogq[g] += qlogq3(q1, q2, q3);
        } else {
            key.clear();
            serialize_eqn(eqn, key);
            auto it = gen_groups.find(key);
            size_t g;
            if (it == gen_groups.end()) {
                g = b.gen_rep.size();
                gen_groups.emplace(key, g);
                b.gen_rep.push_back(qdata[j].get());
                b.gen_q.push_back({0.0, 0.0, 0.0});
                b.gen_qlogq.push_back(0.0);
                std::vector<int> ps;
                eqn_params(eqn, ps);
                b.gen_params.push_back(std::move(ps));
            } else {
                g = it->second;
            }
            b.gen_q[g][0] += q1;
            b.gen_q[g][1] += q2;
            b.gen_q[g][2] += q3;
            b.gen_qlogq[g] += qlogq3(q1, q2, q3);
        }
    }
    b.off.reserve(tree_keys.size() + 1);
    for (const auto& k : tree_keys) {
        b.off.push_back(b.coal.size());
        b.coal.insert(b.coal.end(), k.begin(), k.end());
    }
    b.off.push_back(b.coal.size());
    return b;
}

double computelossandgradient(const TreeQuartetBatch& batch, const std::vector<QuartetData>& /*qdata*/,
                              const double* params, size_t nparam, double* grad, const CFMatrix& /*q*/, double alpha,
                              LossWorkspace& ws, int nthreads) {
    const int np = static_cast<int>(nparam);
    nthreads = std::max(1, nthreads);
    LossWorkspace::Impl* w = ws.get(np, nthreads);
    for (size_t p = 0; p < nparam; ++p) grad[p] = 0.0;
    const double log3 = std::log(3.0);

    // --- tree-like groups ---
    const size_t ng = batch.qres.size();
    std::fill(w->tree_grad.begin(), w->tree_grad.end(), 0.0);
    std::fill(w->partial.begin(), w->partial.end(), 0.0);
    parallel_for(ng, nthreads, [&](size_t lo, size_t hi, int tid) {
        double* g = &w->tree_grad[static_cast<size_t>(tid) * nparam];
        double total = 0.0;
        for (size_t j = lo; j < hi; ++j) {
            const size_t a = batch.off[j], b = batch.off[j + 1];
            double branchsum = 0.0;
            for (size_t k = a; k < b; ++k) branchsum += params[batch.coal[k]];
            const double E = std::exp(-branchsum);
            const double res = 1.0 - (2.0 / 3.0) * E;
            const double oth = E / 3.0;
            const double res_c = std::max(res, 1e-9);
            const double oth_c = std::max(oth, 1e-9);
            const double qr = batch.qres[j], qo = batch.qoth[j];
            const double logoth = (oth >= 1e-9) ? (-branchsum - log3) : std::log(1e-9);
            total += qr * std::log(res_c) + qo * logoth - batch.qlogq[j];
            const double gc = qr * ((2.0 / 3.0) * E) / res_c - qo * ((1.0 / 3.0) * E) / oth_c;
            for (size_t k = a; k < b; ++k) g[batch.coal[k]] += gc;
        }
        w->partial[static_cast<size_t>(tid)] += total;
    });
    double total = 0.0;
    for (int t = 0; t < nthreads; ++t) {
        total += w->partial[static_cast<size_t>(t)];
        const double* g = &w->tree_grad[static_cast<size_t>(t) * nparam];
        for (size_t p = 0; p < nparam; ++p) grad[p] += g[p];
    }

    // --- reticulate groups ---
    if (!batch.gen_rep.empty()) {
        for (int t = 0; t < nthreads; ++t) {
            auto& tg = w->scratch[static_cast<size_t>(t)].total_grad;
            std::fill(tg.begin(), tg.end(), 0.0);
            w->partial[static_cast<size_t>(t)] = 0.0;
        }
        parallel_for(batch.gen_rep.size(), nthreads, [&](size_t lo, size_t hi, int tid) {
            ThreadScratch& ts = w->scratch[static_cast<size_t>(tid)];
            double* ig = ts.quartet_grad.data();
            double tot = 0.0;
            for (size_t g = lo; g < hi; ++g) {
                const auto& ps = batch.gen_params[g];
                for (int p : ps) ig[3 * p] = ig[3 * p + 1] = ig[3 * p + 2] = 0.0;
                auto r = recur(batch.gen_rep[g]->eqn, params, ig, 1.0, alpha);
                double e1 = std::max(r[0], 1e-9), e2 = std::max(r[1], 1e-9);
                double e3 = std::max(1.0 - r[0] - r[1], 1e-9);
                const auto& Q = batch.gen_q[g];
                tot += Q[0] * std::log(e1) + Q[1] * std::log(e2) + Q[2] * std::log(e3) - batch.gen_qlogq[g];
                const double i1 = Q[0] / e1, i2 = Q[1] / e2, i3 = Q[2] / e3;
                for (int p : ps)
                    ts.total_grad[static_cast<size_t>(p)] += ig[3 * p] * i1 + ig[3 * p + 1] * i2 + ig[3 * p + 2] * i3;
            }
            w->partial[static_cast<size_t>(tid)] = tot;
        });
        for (int t = 0; t < nthreads; ++t) {
            total += w->partial[static_cast<size_t>(t)];
            const auto& tg = w->scratch[static_cast<size_t>(t)].total_grad;
            for (size_t p = 0; p < nparam; ++p) grad[p] += tg[p];
        }
    }
    return total;
}

void preparebatch(ClassBatch& batch) {
    batch.dag.build(batch.gen);
    batch.root.clear();
    batch.rootQ.clear();
    std::vector<int32_t> at(batch.dag.node.size(), -1);
    for (size_t g = 0; g < batch.gen.size(); ++g) {
        const uint32_t r = batch.dag.root[g];
        if (at[r] < 0) {
            at[r] = static_cast<int32_t>(batch.root.size());
            batch.root.push_back(r);
            batch.rootQ.push_back(batch.genQ[g]);
        } else {
            auto& Q = batch.rootQ[static_cast<size_t>(at[r])];
            for (size_t k = 0; k < 3; ++k) Q[k] += batch.genQ[g][k];
        }
    }
}

double evalclassvalues(const ClassBatch& batch, const double* params, double alpha, LossWorkspace& ws,
                       int nthreads) {
    nthreads = std::max(1, nthreads);
    LossWorkspace::Impl* w = ws.get(batch.npid, nthreads);
    const double log3 = std::log(3.0);

    // --- tree-like classes ---
    const size_t ng = batch.qres.size();
    w->treeE.resize(ng);
    std::fill(w->partial.begin(), w->partial.end(), 0.0);
    parallel_for(ng, nthreads, [&](size_t lo, size_t hi, int tid) {
        double total = 0.0;
        for (size_t j = lo; j < hi; ++j) {
            const size_t a = batch.off[j], b = batch.off[j + 1];
            double branchsum = 0.0;
            for (size_t k = a; k < b; ++k) branchsum += params[batch.coal[k]];
            const double E = std::exp(-branchsum);
            w->treeE[j] = E;
            const double res = 1.0 - (2.0 / 3.0) * E;
            const double oth = E / 3.0;
            const double logoth = (oth >= 1e-9) ? (-branchsum - log3) : std::log(1e-9);
            total += batch.qres[j] * std::log(std::max(res, 1e-9)) + batch.qoth[j] * logoth;
        }
        w->partial[static_cast<size_t>(tid)] += total;
    });
    double total = -batch.qlogq;
    for (int t = 0; t < nthreads; ++t) total += w->partial[static_cast<size_t>(t)];

    // --- reticulate classes ---
    if (!batch.gen.empty()) {
        if (batch.dag.root.size() != batch.gen.size()) throw SnaqError("evalclasses: preparebatch was not called");
        batch.dag.values(params, alpha, w->dag);
        const double* val = w->dag.val.data();
        for (size_t g = 0; g < batch.root.size(); ++g) {
            const double* v = val + 2 * static_cast<size_t>(batch.root[g]);
            const double e1 = std::max(v[0], 1e-9), e2 = std::max(v[1], 1e-9);
            const double e3 = std::max(1.0 - v[0] - v[1], 1e-9);
            const auto& Q = batch.rootQ[g];
            total += Q[0] * std::log(e1) + Q[1] * std::log(e2) + Q[2] * std::log(e3);
        }
    }
    return total;
}

void evalclassderivatives(const ClassBatch& batch, const double* params, double* grad, double alpha,
                          LossWorkspace& ws, int nthreads, const int* dense, int ndense, double* hess) {
    const size_t nparam = static_cast<size_t>(batch.npid);
    nthreads = hess ? 1 : std::max(1, nthreads);
    LossWorkspace::Impl* w = ws.get(batch.npid, nthreads);
    for (size_t p = 0; p < nparam; ++p) grad[p] = 0.0;
    const size_t nd = static_cast<size_t>(ndense);
    // The curvature is symmetric: each pair of parameters is added once, on
    // one side of the diagonal, and the two sides are summed at the end.
    if (hess) std::fill(hess, hess + nd * nd, 0.0);

    // --- tree-like classes ---
    const size_t ng = batch.qres.size();
    if (w->treeE.size() != ng) throw SnaqError("evalclassderivatives: evalclassvalues was not called");
    std::fill(w->tree_grad.begin(), w->tree_grad.end(), 0.0);
    parallel_for(ng, nthreads, [&](size_t lo, size_t hi, int tid) {
        double* g = &w->tree_grad[static_cast<size_t>(tid) * nparam];
        for (size_t j = lo; j < hi; ++j) {
            const size_t a = batch.off[j], b = batch.off[j + 1];
            const double E = w->treeE[j];
            const double res_c = std::max(1.0 - (2.0 / 3.0) * E, 1e-9);
            const double oth_c = std::max(E / 3.0, 1e-9);
            const double qr = batch.qres[j], qo = batch.qoth[j];
            const double gc = qr * ((2.0 / 3.0) * E) / res_c - qo * ((1.0 / 3.0) * E) / oth_c;
            for (size_t k = a; k < b; ++k) g[batch.coal[k]] += gc;
            if (hess) {
                const double curv = qr * ((2.0 / 3.0) * E) / (res_c * res_c);
                for (size_t k = a; k < b; ++k) {
                    const int dk = dense[batch.coal[k]];
                    if (dk < 0) continue;
                    double* row = hess + static_cast<size_t>(dk) * nd;
                    for (size_t l = a; l <= k; ++l) {
                        const int dl = dense[batch.coal[l]];
                        if (dl >= 0) row[dl] += curv;
                    }
                }
            }
        }
    });
    for (int t = 0; t < nthreads; ++t) {
        const double* g = &w->tree_grad[static_cast<size_t>(t) * nparam];
        for (size_t p = 0; p < nparam; ++p) grad[p] += g[p];
    }

    // --- reticulate classes ---
    if (!batch.gen.empty()) {
        const EqnDag& dag = batch.dag;
        dag.derivatives(params, alpha, w->dag);
        const double* val = w->dag.val.data();
        const double* jac = w->dag.jac.data();
        for (size_t g = 0; g < batch.root.size(); ++g) {
            const EqnDag::Node& n = dag.node[batch.root[g]];
            const double* v = val + 2 * static_cast<size_t>(batch.root[g]);
            const double e1 = std::max(v[0], 1e-9), e2 = std::max(v[1], 1e-9);
            const double e3 = std::max(1.0 - v[0] - v[1], 1e-9);
            const auto& Q = batch.rootQ[g];
            const double i1 = Q[0] / e1, i2 = Q[1] / e2, i3 = Q[2] / e3;
            // d e3 = -(d e1 + d e2)
            const double a1 = i1 - i3, a2 = i2 - i3;
            const double* J = jac + 2 * n.poff;
            const int* ps = dag.pid.data() + n.poff;
            for (uint32_t t = 0; t < n.np; ++t) grad[ps[t]] += a1 * J[2 * t] + a2 * J[2 * t + 1];
            if (hess) {
                // sum over k of (Q_k / e_k^2) d e_k d e_k', as a quadratic form in (d e1, d e2)
                const double c3 = i3 / e3, m11 = i1 / e1 + c3, m22 = i2 / e2 + c3;
                if (w->jidx.size() < n.np) {
                    w->jidx.resize(n.np);
                    for (auto* x : {&w->j1, &w->j2, &w->u1, &w->u2}) x->resize(n.np);
                }
                int* idx = w->jidx.data();
                double *j1 = w->j1.data(), *j2 = w->j2.data(), *u1 = w->u1.data(), *u2 = w->u2.data();
                uint32_t m = 0;
                for (uint32_t t = 0; t < n.np; ++t) {
                    const int d = dense[ps[t]];
                    if (d < 0) continue;
                    idx[m] = d;
                    j1[m] = J[2 * t];
                    j2[m] = J[2 * t + 1];
                    u1[m] = m11 * j1[m] + c3 * j2[m];
                    u2[m] = c3 * j1[m] + m22 * j2[m];
                    ++m;
                }
                for (uint32_t a = 0; a < m; ++a) {
                    double* row = hess + static_cast<size_t>(idx[a]) * nd;
                    const double x1 = j1[a], x2 = j2[a];
                    for (uint32_t b = 0; b <= a; ++b) row[idx[b]] += x1 * u1[b] + x2 * u2[b];
                }
            }
        }
    }
    if (hess) {
        for (size_t a = 0; a < nd; ++a)
            for (size_t b = 0; b < a; ++b) {
                const double sum = hess[a * nd + b] + hess[b * nd + a];
                hess[a * nd + b] = hess[b * nd + a] = sum;
            }
    }
}

double evalclasses(const ClassBatch& batch, const double* params, double* grad, double alpha, LossWorkspace& ws,
                   int nthreads, const int* dense, int ndense, double* hess) {
    const double f = evalclassvalues(batch, params, alpha, ws, nthreads);
    evalclassderivatives(batch, params, grad, alpha, ws, nthreads, dense, ndense, hess);
    return f;
}

double evalclasses_reference(const ClassBatch& batch, const double* params, double* grad, double alpha,
                             LossWorkspace& ws, int nthreads, const int* dense, int ndense, double* hess) {
    const int np = batch.npid;
    const size_t nparam = static_cast<size_t>(np);
    nthreads = hess ? 1 : std::max(1, nthreads);
    LossWorkspace::Impl* w = ws.get(np, nthreads);
    for (size_t p = 0; p < nparam; ++p) grad[p] = 0.0;
    if (hess) std::fill(hess, hess + static_cast<size_t>(ndense) * static_cast<size_t>(ndense), 0.0);
    const double log3 = std::log(3.0);

    // --- tree-like classes ---
    const size_t ng = batch.qres.size();
    std::fill(w->tree_grad.begin(), w->tree_grad.end(), 0.0);
    std::fill(w->partial.begin(), w->partial.end(), 0.0);
    parallel_for(ng, nthreads, [&](size_t lo, size_t hi, int tid) {
        double* g = &w->tree_grad[static_cast<size_t>(tid) * nparam];
        double total = 0.0;
        for (size_t j = lo; j < hi; ++j) {
            const size_t a = batch.off[j], b = batch.off[j + 1];
            double branchsum = 0.0;
            for (size_t k = a; k < b; ++k) branchsum += params[batch.coal[k]];
            const double E = std::exp(-branchsum);
            const double res = 1.0 - (2.0 / 3.0) * E;
            const double oth = E / 3.0;
            const double res_c = std::max(res, 1e-9);
            const double oth_c = std::max(oth, 1e-9);
            const double qr = batch.qres[j], qo = batch.qoth[j];
            const double logoth = (oth >= 1e-9) ? (-branchsum - log3) : std::log(1e-9);
            total += qr * std::log(res_c) + qo * logoth;
            const double gc = qr * ((2.0 / 3.0) * E) / res_c - qo * ((1.0 / 3.0) * E) / oth_c;
            for (size_t k = a; k < b; ++k) g[batch.coal[k]] += gc;
            if (hess) {
                const double curv = qr * ((2.0 / 3.0) * E) / (res_c * res_c);
                for (size_t k = a; k < b; ++k) {
                    const int dk = dense[batch.coal[k]];
                    if (dk < 0) continue;
                    double* row = hess + static_cast<size_t>(dk) * static_cast<size_t>(ndense);
                    for (size_t l = a; l < b; ++l) {
                        const int dl = dense[batch.coal[l]];
                        if (dl >= 0) row[dl] += curv;
                    }
                }
            }
        }
        w->partial[static_cast<size_t>(tid)] += total;
    });
    double total = -batch.qlogq;
    for (int t = 0; t < nthreads; ++t) {
        total += w->partial[static_cast<size_t>(t)];
        const double* g = &w->tree_grad[static_cast<size_t>(t) * nparam];
        for (size_t p = 0; p < nparam; ++p) grad[p] += g[p];
    }

    // --- reticulate classes ---
    if (!batch.gen.empty()) {
        for (int t = 0; t < nthreads; ++t) {
            auto& tg = w->scratch[static_cast<size_t>(t)].total_grad;
            std::fill(tg.begin(), tg.end(), 0.0);
            w->partial[static_cast<size_t>(t)] = 0.0;
        }
        parallel_for(batch.gen.size(), nthreads, [&](size_t lo, size_t hi, int tid) {
            ThreadScratch& ts = w->scratch[static_cast<size_t>(tid)];
            double* ig = ts.quartet_grad.data();
            double tot = 0.0;
            std::vector<int> ps;  // the parameters of the class
            for (size_t g = lo; g < hi; ++g) {
                ps.clear();
                eqn_params(batch.gen[g]->eqn, ps);
                for (int p : ps) ig[3 * p] = ig[3 * p + 1] = ig[3 * p + 2] = 0.0;
                auto r = recur(batch.gen[g]->eqn, params, ig, 1.0, alpha);
                double e1 = std::max(r[0], 1e-9), e2 = std::max(r[1], 1e-9);
                double e3 = std::max(1.0 - r[0] - r[1], 1e-9);
                const auto& Q = batch.genQ[g];
                tot += Q[0] * std::log(e1) + Q[1] * std::log(e2) + Q[2] * std::log(e3);
                const double i1 = Q[0] / e1, i2 = Q[1] / e2, i3 = Q[2] / e3;
                for (int p : ps)
                    ts.total_grad[static_cast<size_t>(p)] += ig[3 * p] * i1 + ig[3 * p + 1] * i2 + ig[3 * p + 2] * i3;
                if (hess) {
                    const double c1 = i1 / e1, c2 = i2 / e2, c3 = i3 / e3;
                    for (int p : ps) {
                        const int dp = dense[p];
                        if (dp < 0) continue;
                        const double a1 = c1 * ig[3 * p], a2 = c2 * ig[3 * p + 1], a3 = c3 * ig[3 * p + 2];
                        double* row = hess + static_cast<size_t>(dp) * static_cast<size_t>(ndense);
                        for (int q : ps) {
                            const int dq = dense[q];
                            if (dq >= 0) row[dq] += a1 * ig[3 * q] + a2 * ig[3 * q + 1] + a3 * ig[3 * q + 2];
                        }
                    }
                }
            }
            w->partial[static_cast<size_t>(tid)] = tot;
        });
        for (int t = 0; t < nthreads; ++t) {
            total += w->partial[static_cast<size_t>(t)];
            const auto& tg = w->scratch[static_cast<size_t>(t)].total_grad;
            for (size_t p = 0; p < nparam; ++p) grad[p] += tg[p];
        }
    }
    return total;
}

std::array<double, 3> classexpectedCF(const ClassEqn& ce, const double* params, int /*npid*/, double alpha) {
    const auto r = recur(ce.eqn, params, nullptr, 1.0, alpha);
    return {r[0], r[1], 1.0 - r[0] - r[1]};
}

double treeclassbound(double qr, double qo, bool hasparams) {
    static const double log3 = std::log(3.0), logfloor = std::log(1e-9);
    if (!hasparams) return -(qr + qo) * log3;  // the path length is 0
    // maximize qr log(1 - 2z/3) + qo log(z/3) over z = exp(-s) in (0, 1]
    if (qo <= 0.0) return 0.0;
    const double z = 3.0 * qo / (2.0 * (qr + qo));
    const double best = z >= 1.0 ? -(qr + qo) * log3 : qr * std::log(1.0 - 2.0 * z / 3.0) + qo * std::log(z / 3.0);
    return std::max(best, qo * logfloor);  // the loss floors expected CFs at 1e-9
}

double generalclassbound(const std::array<double, 3>& Q) {
    static const double logfloor = std::log(1e-9);
    const double n = Q[0] + Q[1] + Q[2];
    double ub = 4e-9 * n;  // expected CFs floored at 1e-9 can sum to slightly more than 1
    for (int k = 0; k < 3; ++k)
        if (Q[k] > 0.0) ub += Q[k] * std::max(std::log(Q[k] / n), logfloor);
    return ub;
}

double computelossandgradient(const std::vector<QuartetData>& qdata, const std::vector<double>& params,
                              std::vector<double>& grad, const CFMatrix& q, double alpha, int nthreads) {
    const int np = static_cast<int>(params.size());
    nthreads = std::max(1, nthreads);
    grad.assign(params.size(), 0.0);
    std::vector<ThreadScratch> scratch;
    for (int t = 0; t < nthreads; ++t) scratch.emplace_back(np);
    std::vector<double> partial(static_cast<size_t>(nthreads), 0.0);
    parallel_for(qdata.size(), nthreads, [&](size_t lo, size_t hi, int tid) {
        ThreadScratch& ts = scratch[static_cast<size_t>(tid)];
        double tot = 0.0;
        for (size_t j = lo; j < hi; ++j) {
            const QuartetDataImpl& qd = *qdata[j];
            if (qd.eqn.division_H == -1)
                tot += treequartetlossandgradient(qd.eqn, params.data(), q.row(j), ts.total_grad.data());
            else
                tot += quartetlossandgradient_general(qd, params.data(), q.row(j), alpha, ts);
        }
        partial[static_cast<size_t>(tid)] = tot;
    });
    double total = 0.0;
    for (int t = 0; t < nthreads; ++t) {
        total += partial[static_cast<size_t>(t)];
        for (size_t p = 0; p < params.size(); ++p) grad[p] += scratch[static_cast<size_t>(t)].total_grad[p];
    }
    return total;
}

double computeSNaQscore(const std::vector<QuartetData>& qdata, const std::vector<double>& params, const CFMatrix& q,
                        double rho, int nthreads) {
    std::vector<double> grad;
    return computelossandgradient(qdata, params, grad, q, rhotoalpha(rho), nthreads);
}

std::array<double, 2> computeexpectedCF(const QuartetDataImpl& qd, const std::vector<double>& params, double rho) {
    const auto r = recur(qd.eqn, params.data(), nullptr, 1.0, rhotoalpha(rho));
    return {r[0], r[1]};
}

}  // namespace snaqpp

namespace snaqpp {

CFMatrix computeexpectedCFmatrix(Network& net, double rho, int nthreads) {
    fixnegativeedges(net);
    checkgammas(net);
    const size_t nq = static_cast<size_t>(nchoose4(net.numtaxa));
    CFMatrix dummy;
    dummy.v.assign(3 * nq, 1.0 / 3.0);
    ClassOptions co;
    co.perquartet = true;
    co.nthreads = nthreads;
    QuartetClasses cls = buildclasses(net, QuartetList::allquartets(net.numtaxa), dummy, co);
    const std::vector<double> params = pidparams(net, cls.npid);
    const double alpha = rhotoalpha(rho);
    std::vector<std::array<double, 3>> ecf(cls.size());
    for (size_t c = 0; c < cls.size(); ++c) ecf[c] = classexpectedCF(*cls.eqn[c], params.data(), cls.npid, alpha);
    CFMatrix m;
    m.v.resize(3 * nq);
    for (size_t j = 0; j < nq; ++j) {
        const auto& e = ecf[cls.qclass[j]];
        const uint8_t cols = cls.qcols[j];
        for (int k = 0; k < 3; ++k) m(j, (cols >> (2 * k)) & 3) = e[static_cast<size_t>(k)];
    }
    return m;
}

double computeSNaQscore_net(Network& N, const CFMatrix& q, double rho, int nthreads) {
    Network c = N.deepcopynetwork();
    semidirectnetwork(c);
    fixnegativeedges(c);
    checkgammas(c);
    ClassOptions co;
    co.qlogq = obsqlogq(q);
    co.nthreads = nthreads;
    QuartetClasses cls = buildclasses(c, QuartetList::allquartets(c.numtaxa), q, co);
    double s = scoreclasses(c, cls, rho, nthreads);
    N.fscore = s;
    return s;
}

}  // namespace snaqpp
