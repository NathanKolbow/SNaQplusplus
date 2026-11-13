#include "likelihood.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <unordered_map>

#include "classes.hpp"
#include "parallel.hpp"
#include "snaq_network.hpp"

namespace csnaq {

double rhotoalpha(double rho) { return rho == 0.0 ? std::numeric_limits<double>::infinity() : (1.0 - rho) / rho; }
double alphatorho(double alpha) { return std::isinf(alpha) ? 0.0 : 1.0 / (1.0 + alpha); }

namespace {

// ------------------------------------------------------------------------
// RunningGradient (running_gradient.jl)
// ------------------------------------------------------------------------
struct RunningGradient {
    double shared[3] = {1.0, 1.0, 1.0};
    std::vector<int> exceptionparam;
    std::vector<double> exceptionvalue;  // 3 x capacity, column-major (3*s + k)
    std::vector<int> exceptionslot;      // nparam: 0 if not an exception, else 1-based slot
    std::vector<char> params_seen;
    int nexceptions = 0;
    std::vector<double> savestack;
    size_t savetop = 0;

    explicit RunningGradient(int nparam)
        : exceptionparam(8, 0),
          exceptionvalue(3 * 8, 0.0),
          exceptionslot(static_cast<size_t>(nparam), 0),
          params_seen(static_cast<size_t>(nparam), 0),
          savestack(64, 0.0) {}

    void reset() {
        shared[0] = shared[1] = shared[2] = 1.0;
        for (int s = 0; s < nexceptions; ++s) {
            exceptionslot[static_cast<size_t>(exceptionparam[static_cast<size_t>(s)])] = 0;
            params_seen[static_cast<size_t>(exceptionparam[static_cast<size_t>(s)])] = 0;
        }
        nexceptions = 0;
        savetop = 0;
    }
    double value(int p, int k) const {
        int s = exceptionslot[static_cast<size_t>(p)];
        return s == 0 ? shared[k] : exceptionvalue[3 * static_cast<size_t>(s - 1) + static_cast<size_t>(k)];
    }
    double& exval(int s0, int k) { return exceptionvalue[3 * static_cast<size_t>(s0) + static_cast<size_t>(k)]; }
    void addexception(int p) {
        if (exceptionslot[static_cast<size_t>(p)] != 0) return;
        if (nexceptions == static_cast<int>(exceptionparam.size())) {
            exceptionparam.resize(exceptionparam.size() * 2);
            exceptionvalue.resize(exceptionvalue.size() * 2);
        }
        int s0 = nexceptions++;
        exceptionparam[static_cast<size_t>(s0)] = p;
        exceptionslot[static_cast<size_t>(p)] = s0 + 1;
        exval(s0, 0) = shared[0];
        exval(s0, 1) = shared[1];
        exval(s0, 2) = shared[2];
    }
    void dropexceptionsto(int nex0) {
        for (int s = nexceptions - 1; s >= nex0; --s)
            exceptionslot[static_cast<size_t>(exceptionparam[static_cast<size_t>(s)])] = 0;
        nexceptions = nex0;
    }
    size_t save() {
        size_t need = savetop + 3 + 3 * static_cast<size_t>(nexceptions);
        if (savestack.size() < need) savestack.resize(std::max(need, 2 * savestack.size()));
        savestack[savetop] = shared[0];
        savestack[savetop + 1] = shared[1];
        savestack[savetop + 2] = shared[2];
        size_t o = savetop + 3;
        for (int s = 0; s < nexceptions; ++s) {
            savestack[o] = exval(s, 0);
            savestack[o + 1] = exval(s, 1);
            savestack[o + 2] = exval(s, 2);
            o += 3;
        }
        savetop = need;
        return need;
    }
    void restore(size_t top) {
        size_t base = top - 3 - 3 * static_cast<size_t>(nexceptions);
        shared[0] = savestack[base];
        shared[1] = savestack[base + 1];
        shared[2] = savestack[base + 2];
        size_t o = base + 3;
        for (int s = 0; s < nexceptions; ++s) {
            exval(s, 0) = savestack[o];
            exval(s, 1) = savestack[o + 1];
            exval(s, 2) = savestack[o + 2];
            o += 3;
        }
        savetop = base;
    }
    void scale(double f, int skip = -1) {
        shared[0] *= f;
        shared[1] *= f;
        shared[2] *= f;
        for (int s = 0; s < nexceptions; ++s) {
            if (exceptionparam[static_cast<size_t>(s)] == skip) continue;
            exval(s, 0) *= f;
            exval(s, 1) *= f;
            exval(s, 2) *= f;
        }
    }
    void scaleexception(int p, double f) {
        int s = exceptionslot[static_cast<size_t>(p)];
        if (s == 0) return;
        exval(s - 1, 0) *= f;
        exval(s - 1, 1) *= f;
        exval(s - 1, 2) *= f;
    }
};

inline double coalbranchsum(const RecursiveCFEquation& eqn, const double* params) {
    double s = 0.0;
    for (int p : eqn.coal_edges) s += params[p];
    return s;
}

template <class V>
inline bool iscoaledge(const V& coal, int p) {
    for (int e : coal)
        if (e == p) return true;
    return false;
}

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

// gs: nparam x 3, row-major (gs[3*p + k])
std::array<double, 2> recur(const RecursiveCFEquation& eqn, const double* params, double* gs, RunningGradient& rg,
                            double alpha) {
    const auto& coal = eqn.coal_edges;
    const int wc = eqn.which_coal;
    if (eqn.division_H == -1) {
        const double E = std::exp(-coalbranchsum(eqn, params));
        const double d[3] = {wc == 1 ? 2.0 / 3.0 * E : -1.0 / 3.0 * E, wc == 2 ? 2.0 / 3.0 * E : -1.0 / 3.0 * E,
                             wc == 3 ? 2.0 / 3.0 * E : -1.0 / 3.0 * E};
        const double f[3] = {wc == 1 ? 1 - 2.0 / 3.0 * E : 1.0 / 3.0 * E, wc == 2 ? 1 - 2.0 / 3.0 * E : 1.0 / 3.0 * E,
                             wc == 3 ? 1 - 2.0 / 3.0 * E : 1.0 / 3.0 * E};
        for (int p : coal)
            for (int k = 0; k < 3; ++k) gs[3 * p + k] += rg.value(p, k) * d[k];
        for (int s = 0; s < rg.nexceptions; ++s) {
            int p = rg.exceptionparam[static_cast<size_t>(s)];
            if (!rg.params_seen[static_cast<size_t>(p)]) continue;
            if (iscoaledge(coal, p)) continue;
            for (int k = 0; k < 3; ++k) gs[3 * p + k] += rg.exval(s, k) * f[k];
        }
        if (wc == 1) return {1 - 2.0 / 3.0 * E, 1.0 / 3.0 * E};
        if (wc == 2) return {1.0 / 3.0 * E, 1 - 2.0 / 3.0 * E};
        return {1.0 / 3.0 * E, 1.0 / 3.0 * E};
    } else if (eqn.divisions.size() == 4) {
        double e1 = 0.0, e2 = 0.0;
        const int nex0 = rg.nexceptions;
        double early = std::exp(-coalbranchsum(eqn, params));
        early = std::max(early, 1e-9);
        const int wk = wc - 1;
        if (eqn.can_coalesce_here) {
            if (wc == 1)
                e1 += 1 - early;
            else if (wc == 2)
                e2 += 1 - early;
            for (int p : coal) {
                if (rg.params_seen[static_cast<size_t>(p)]) throw SnaqError("Already seen this param??");
                gs[3 * p + wk] += rg.value(p, wk) * early;
            }
            for (int s = 0; s < rg.nexceptions; ++s) {
                int p = rg.exceptionparam[static_cast<size_t>(s)];
                if (!rg.params_seen[static_cast<size_t>(p)]) continue;
                gs[3 * p + wk] += rg.exval(s, wk) * (1 - early);
            }
            for (int p : coal) rg.params_seen[static_cast<size_t>(p)] = 1;
        }
        for (int p : coal) rg.addexception(p);
        rg.scale(early);
        const int H = eqn.division_H;
        if (rg.params_seen[static_cast<size_t>(H)]) throw SnaqError("Already seen this param??");
        rg.addexception(H);
        rg.params_seen[static_cast<size_t>(H)] = 1;
        const double g = params[H];
        for (int di = 1; di <= 4; ++di) {
            double split_grad = quadsplitprobabilitygradient(di, g, alpha);
            double split_prob = quadsplitprobability(di, g, alpha);
            size_t top = rg.save();
            rg.scale(split_prob, H);
            rg.scaleexception(H, split_grad);
            for (int e : coal) rg.scaleexception(e, -1.0);
            auto r = recur(eqn.divisions[static_cast<size_t>(di - 1)], params, gs, rg, alpha);
            rg.restore(top);
            double fct = early * split_prob;
            e1 += fct * r[0];
            e2 += fct * r[1];
        }
        rg.scale(1.0 / early);
        for (int e : coal) rg.params_seen[static_cast<size_t>(e)] = 0;
        rg.params_seen[static_cast<size_t>(H)] = 0;
        rg.dropexceptionsto(nex0);
        return {e1, e2};
    } else {
        if (eqn.can_coalesce_here) throw SnaqError("Can coalesce w/ 1 taxa splitting at a single hybrid??");
        const int nex0 = rg.nexceptions;
        const int H = eqn.division_H;
        double g = params[H];
        g = std::max(1e-9, g);
        g = std::min(1.0 - 1e-9, g);
        rg.addexception(H);
        rg.params_seen[static_cast<size_t>(H)] = 1;

        rg.scale(g, H);
        auto r = recur(eqn.divisions[0], params, gs, rg, alpha);
        double e1 = g * r[0], e2 = g * r[1];
        rg.scale(1.0 / g, H);

        rg.scale(1 - g, H);
        rg.scaleexception(H, -1.0);
        auto s = recur(eqn.divisions[1], params, gs, rg, alpha);
        rg.scale(1.0 / (1 - g), H);
        rg.scaleexception(H, -1.0);

        rg.params_seen[static_cast<size_t>(H)] = 0;
        rg.dropexceptionsto(nex0);
        e1 += (1 - g) * s[0];
        e2 += (1 - g) * s[1];
        return {e1, e2};
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
    RunningGradient rg;
    std::vector<double> total_grad;  // nparam
    explicit ThreadScratch(int np)
        : quartet_grad(3 * static_cast<size_t>(np), 0.0), rg(np), total_grad(static_cast<size_t>(np), 0.0) {}
};

// General quartet: loss, gradient accumulated into ts.total_grad.
double quartetlossandgradient_general(const QuartetDataImpl& qd, const double* params, const double* q, double alpha,
                                      ThreadScratch& ts) {
    double* ig = ts.quartet_grad.data();
    ts.rg.reset();
    for (int p : qd.relevant_params) {
        ig[3 * p] = 0.0;
        ig[3 * p + 1] = 0.0;
        ig[3 * p + 2] = 0.0;
    }
    auto r = recur(qd.eqn, params, ig, ts.rg, alpha);
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
                ts.rg.reset();
                for (int p : ps) ig[3 * p] = ig[3 * p + 1] = ig[3 * p + 2] = 0.0;
                auto r = recur(batch.gen_rep[g]->eqn, params, ig, ts.rg, alpha);
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

double evalclasses(const ClassBatch& batch, const double* params, double* grad, double alpha, LossWorkspace& ws,
                   int nthreads, const int* dense, int ndense, double* hess) {
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
            for (size_t g = lo; g < hi; ++g) {
                const auto& ps = batch.gen[g]->params;
                ts.rg.reset();
                for (int p : ps) ig[3 * p] = ig[3 * p + 1] = ig[3 * p + 2] = 0.0;
                auto r = recur(batch.gen[g]->eqn, params, ig, ts.rg, alpha);
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

std::array<double, 3> classexpectedCF(const ClassEqn& ce, const double* params, int npid, double alpha) {
    static thread_local std::vector<double> gs;
    static thread_local std::unique_ptr<RunningGradient> rg;
    static thread_local int rgsize = -1;
    if (gs.size() < 3 * static_cast<size_t>(npid) + 3) gs.assign(3 * static_cast<size_t>(npid) + 3, 0.0);
    if (rgsize < npid) {
        rg = std::make_unique<RunningGradient>(npid);
        rgsize = npid;
    }
    rg->reset();
    auto r = recur(ce.eqn, params, gs.data(), *rg, alpha);
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
    double alpha = rhotoalpha(rho);
    std::vector<double> gs(3 * params.size() + 3, 0.0);
    RunningGradient rg(static_cast<int>(params.size()));
    return recur(qd.eqn, params.data(), gs.data(), rg, alpha);
}

}  // namespace csnaq

namespace csnaq {

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

}  // namespace csnaq
