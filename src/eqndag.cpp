#include "eqndag.hpp"

#include <algorithm>
#include <cmath>

#include "likelihood.hpp"

namespace csnaq {

namespace {

inline uint64_t mixa(uint64_t x) {
    x ^= x >> 33;
    x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33;
    x *= 0xc4ceb9fe1a85ec53ULL;
    return x ^ (x >> 33);
}
inline uint64_t mixb(uint64_t x) {
    x ^= x >> 30;
    x *= 0xbf58476d1ce4e5b9ULL;
    x ^= x >> 27;
    x *= 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}
struct Hasher {
    uint64_t a = 0x243f6a8885a308d3ULL, b = 0x13198a2e03707344ULL;
    void add(uint64_t v) {
        a = mixa(a ^ (v + 0x9e3779b97f4a7c15ULL));
        b = mixb(b + (v ^ 0xa4093822299f31d0ULL));
    }
};

struct Key {
    uint64_t a, b;
};

// Open-addressing hash table from a 128-bit key to an index.
class KeyTable {
public:
    explicit KeyTable(size_t expected) {
        size_t n = 64;
        while (n < 2 * expected) n *= 2;
        slots_.assign(n, Slot{0, 0, EMPTY});
    }
    static constexpr uint32_t EMPTY = 0xffffffffu;
    // The value of `k`, or EMPTY.
    uint32_t find(const Key& k) const {
        const size_t mask = slots_.size() - 1;
        for (size_t i = static_cast<size_t>(k.a) & mask;; i = (i + 1) & mask) {
            const Slot& s = slots_[i];
            if (s.value == EMPTY) return EMPTY;
            if (s.a == k.a && s.b == k.b) return s.value;
        }
    }
    void insert(const Key& k, uint32_t value) {
        if (2 * (count_ + 1) > slots_.size()) grow();
        place(k, value);
        ++count_;
    }

private:
    struct Slot {
        uint64_t a, b;
        uint32_t value;
    };
    void place(const Key& k, uint32_t value) {
        const size_t mask = slots_.size() - 1;
        size_t i = static_cast<size_t>(k.a) & mask;
        while (slots_[i].value != EMPTY) i = (i + 1) & mask;
        slots_[i] = Slot{k.a, k.b, value};
    }
    void grow() {
        std::vector<Slot> old(2 * slots_.size(), Slot{0, 0, EMPTY});
        old.swap(slots_);
        for (const Slot& s : old)
            if (s.value != EMPTY) place(Key{s.a, s.b}, s.value);
    }
    std::vector<Slot> slots_;
    size_t count_ = 0;
};

using Sorted = SmallVec<int, 24>;
inline void sortededges(const ParamList& coal, Sorted& out) {
    out.clear();
    for (int p : coal) out.push_back(p);
    std::sort(out.begin(), out.end());
}

EqnHash hashnode(const RecursiveCFEquation& e, std::vector<EqnHash>& out) {
    const size_t at = out.size();
    out.emplace_back();
    Hasher h;
    Sorted coal;
    if (e.division_H == -1) {
        sortededges(e.coal_edges, coal);
        h.add(1);
        h.add(static_cast<uint64_t>(e.which_coal));
        for (int p : coal) h.add(static_cast<uint64_t>(p) + 16);
    } else if (e.divisions.size() == 4) {
        sortededges(e.coal_edges, coal);
        const bool here = e.can_coalesce_here && e.which_coal >= 1;
        h.add(4);
        h.add(static_cast<uint64_t>(e.division_H));
        h.add(here ? static_cast<uint64_t>(e.which_coal) : 0);
        h.add(coal.size());
        for (int p : coal) h.add(static_cast<uint64_t>(p) + 16);
        EqnHash c[4];
        for (size_t k = 0; k < 4; ++k) c[k] = hashnode(e.divisions[k], out);
        if (c[3].a < c[2].a || (c[3].a == c[2].a && c[3].b < c[2].b)) std::swap(c[2], c[3]);
        for (size_t k = 0; k < 4; ++k) h.add(c[k].a), h.add(c[k].b);
    } else {
        if (e.divisions.size() != 2) throw SnaqError("equation with neither 2 nor 4 divisions");
        h.add(2);
        h.add(static_cast<uint64_t>(e.division_H));
        for (size_t k = 0; k < 2; ++k) {
            const EqnHash c = hashnode(e.divisions[k], out);
            h.add(c.a), h.add(c.b);
        }
    }
    EqnHash r;
    r.a = h.a;
    r.b = h.b;
    r.size = static_cast<uint32_t>(out.size() - at);
    out[at] = r;
    return r;
}

// Probabilities of the four ways two lineages pass a hybrid node with minor γ
// `g` (both minor, both major, one each, one each), and their derivatives in
// γ: quadsplitprobability and quadsplitprobabilitygradient for the four types.
inline void splitprobabilities(double g, double a, double* sp, double* sg) {
    if (std::isinf(a)) {
        sp[0] = g * g;
        sp[1] = (1 - g) * (1 - g);
        sp[2] = sp[3] = g * (1 - g);
        sg[0] = 2 * g;
        sg[1] = -2 * (1 - g);
        sg[2] = sg[3] = 1 - 2 * g;
    } else if (a == 0.0) {
        sp[0] = g;
        sp[1] = 1 - g;
        sp[2] = sp[3] = 0;
        sg[0] = 1;
        sg[1] = -1;
        sg[2] = sg[3] = 0.0;
    } else {
        sp[0] = g * (1 / (a + 1) + a / (a + 1) * g);
        sp[1] = (1 - g) * (1 / (a + 1) + a / (a + 1) * (1 - g));
        sp[2] = sp[3] = g * (a / (a + 1)) * (1 - g);
        sg[0] = 1 / (a + 1) + 2 * g * a / (a + 1);
        sg[1] = -1 / (a + 1) - 2 * (1 - g) * a / (a + 1);
        sg[2] = sg[3] = (a / (a + 1)) * (1 - 2 * g);
    }
}

struct Builder {
    EqnDag& dag;
    KeyTable ids, sets;
    std::vector<int> tmp;
    // by pid: whether it is among the parameters of the node being added, and where
    std::vector<uint32_t> seen, pos;
    uint32_t stamp = 0;

    int32_t setof(const ParamList& coal) {
        if (coal.empty()) return -1;
        Sorted s;
        sortededges(coal, s);
        Hasher h;
        for (int p : s) h.add(static_cast<uint64_t>(p));
        const Key k{h.a, h.b};
        const uint32_t found = sets.find(k);
        if (found != KeyTable::EMPTY) return static_cast<int32_t>(found);
        const uint32_t id = static_cast<uint32_t>(dag.setoff.size() - 1);
        for (int p : s) dag.setpid.push_back(p);
        dag.setoff.push_back(static_cast<uint32_t>(dag.setpid.size()));
        sets.insert(k, id);
        return static_cast<int32_t>(id);
    }
    // the parameters of node `i`
    const int* paramsof(uint32_t i, uint32_t& n) const {
        const EqnDag::Node& c = dag.node[i];
        n = c.np;
        if (c.kind != 0) return dag.pid.data() + c.poff;
        return c.set < 0 ? nullptr : dag.setpid.data() + dag.setoff[static_cast<size_t>(c.set)];
    }

    uint32_t add(const RecursiveCFEquation& e, const EqnHash* h) {
        const Key key{h->a, h->b};
        const uint32_t found = ids.find(key);
        if (found != KeyTable::EMPTY) return found;
        EqnDag::Node n;
        if (e.division_H == -1) {
            n.kind = 0;
            n.wc = static_cast<uint8_t>(e.which_coal);
            n.set = setof(e.coal_edges);
            n.np = static_cast<uint32_t>(e.coal_edges.size());
        } else {
            const size_t nd = e.divisions.size();
            n.kind = static_cast<uint8_t>(nd);
            n.H = e.division_H;
            if (nd == 4) {
                n.wc = e.can_coalesce_here && e.which_coal >= 1 ? static_cast<uint8_t>(e.which_coal) : 0;
                n.set = setof(e.coal_edges);
            }
            const EqnHash* hc = h + 1;
            for (size_t k = 0; k < nd; ++k) {
                n.child[k] = add(e.divisions[k], hc);
                hc += hc->size;
            }
            // its parameters: γ, its edges, and the parameters of its children
            tmp.clear();
            ++stamp;
            auto collect = [&](const int* p, uint32_t count) {
                for (uint32_t t = 0; t < count; ++t) {
                    const size_t q = static_cast<size_t>(p[t]);
                    if (q >= seen.size()) {
                        seen.resize(2 * q + 16, 0);
                        pos.resize(2 * q + 16, 0);
                    }
                    if (seen[q] == stamp) continue;
                    seen[q] = stamp;
                    tmp.push_back(p[t]);
                }
            };
            collect(&n.H, 1);
            uint32_t ns = 0;
            if (n.set >= 0) {
                ns = dag.setoff[static_cast<size_t>(n.set) + 1] - dag.setoff[static_cast<size_t>(n.set)];
                collect(dag.setpid.data() + dag.setoff[static_cast<size_t>(n.set)], ns);
            }
            for (size_t k = 0; k < nd; ++k) {
                uint32_t nc;
                const int* cp = paramsof(n.child[k], nc);
                collect(cp, nc);
            }
            std::sort(tmp.begin(), tmp.end());
            n.poff = static_cast<uint32_t>(dag.pid.size());
            n.np = static_cast<uint32_t>(tmp.size());
            for (uint32_t t = 0; t < n.np; ++t) {
                pos[static_cast<size_t>(tmp[t])] = t;
                dag.pid.push_back(tmp[t]);
            }
            n.hpos = pos[static_cast<size_t>(n.H)];
            n.moff = static_cast<uint32_t>(dag.map.size());
            for (size_t k = 0; k < nd; ++k) {
                uint32_t nc;
                const int* cp = paramsof(n.child[k], nc);
                for (uint32_t t = 0; t < nc; ++t) dag.map.push_back(pos[static_cast<size_t>(cp[t])]);
            }
            const int* sp = n.set >= 0 ? dag.setpid.data() + dag.setoff[static_cast<size_t>(n.set)] : nullptr;
            for (uint32_t t = 0; t < ns; ++t) dag.map.push_back(pos[static_cast<size_t>(sp[t])]);
        }
        const uint32_t id = static_cast<uint32_t>(dag.node.size());
        dag.node.push_back(n);
        ids.insert(key, id);
        return id;
    }
};

}  // namespace

std::vector<EqnHash> hashequation(const RecursiveCFEquation& eqn) {
    std::vector<EqnHash> out;
    hashnode(eqn, out);
    return out;
}

void EqnDag::build(const std::vector<const ClassEqn*>& classes) {
    node.clear();
    setoff.assign(1, 0);
    setpid.clear();
    pid.clear();
    map.clear();
    root.resize(classes.size());
    if (classes.empty()) return;
    const size_t nc = classes.size();
    node.reserve(2 * nc);
    pid.reserve(24 * nc);
    map.reserve(32 * nc);
    Builder b{*this, KeyTable(2 * nc), KeyTable(nc / 2), {}, {}, {}, 0};
    std::vector<EqnHash> local;
    for (size_t c = 0; c < classes.size(); ++c) {
        const ClassEqn& ce = *classes[c];
        const EqnHash* h = ce.hash.data();
        if (ce.hash.empty()) {
            local = hashequation(ce.eqn);
            h = local.data();
        }
        root[c] = b.add(ce.eqn, h);
    }
}

void EqnDag::values(const double* params, double alpha, Work& w) const {
    const size_t ns = setoff.size() - 1, nn = node.size();
    w.ex.resize(ns);
    w.val.resize(2 * nn);
    for (size_t s = 0; s < ns; ++s) {
        double sum = 0.0;
        for (uint32_t k = setoff[s]; k < setoff[s + 1]; ++k) sum += params[setpid[k]];
        w.ex[s] = std::exp(-sum);
    }
    double* val = w.val.data();
    double sp[4], sg[4];
    for (size_t i = 0; i < nn; ++i) {
        const Node& n = node[i];
        double* v = val + 2 * i;
        if (n.kind == 0) {
            const double E = n.set < 0 ? 1.0 : w.ex[static_cast<size_t>(n.set)];
            const double res = 1 - 2.0 / 3.0 * E, oth = 1.0 / 3.0 * E;
            v[0] = n.wc == 1 ? res : oth;
            v[1] = n.wc == 2 ? res : oth;
        } else if (n.kind == 2) {
            const double g = std::min(1.0 - 1e-9, std::max(1e-9, params[n.H]));
            const double *a = val + 2 * n.child[0], *b = val + 2 * n.child[1];
            v[0] = g * a[0] + (1 - g) * b[0];
            v[1] = g * a[1] + (1 - g) * b[1];
        } else {
            const double early = std::max(n.set < 0 ? 1.0 : w.ex[static_cast<size_t>(n.set)], 1e-9);
            splitprobabilities(params[n.H], alpha, sp, sg);
            double e0 = 0.0, e1 = 0.0;
            for (int k = 0; k < 4; ++k) {
                const double f = early * sp[k];
                const double* c = val + 2 * n.child[k];
                e0 += f * c[0];
                e1 += f * c[1];
            }
            if (n.wc == 1) e0 += 1 - early;
            if (n.wc == 2) e1 += 1 - early;
            v[0] = e0;
            v[1] = e1;
        }
    }
}

void EqnDag::derivatives(const double* params, double alpha, Work& w) const {
    const size_t nn = node.size();
    w.jac.resize(2 * pid.size());
    const double* val = w.val.data();
    double* jac = w.jac.data();
    double sp[4], sg[4];
    for (size_t i = 0; i < nn; ++i) {
        const Node& n = node[i];
        if (n.kind == 0) continue;
        double* J = jac + 2 * n.poff;
        std::fill(J, J + 2 * n.np, 0.0);
        const uint32_t* m = map.data() + n.moff;
        // adds f * (derivatives of child c)
        auto addchild = [&](uint32_t c, double f) {
            const Node& ch = node[c];
            if (ch.kind == 0) {
                if (ch.set < 0) return;
                const double oth = f * (1.0 / 3.0 * w.ex[static_cast<size_t>(ch.set)]);
                const double d0 = ch.wc == 1 ? 2.0 * oth : -oth, d1 = ch.wc == 2 ? 2.0 * oth : -oth;
                for (uint32_t t = 0; t < ch.np; ++t) {
                    J[2 * m[t]] += d0;
                    J[2 * m[t] + 1] += d1;
                }
            } else {
                const double* Jc = jac + 2 * ch.poff;
                for (uint32_t t = 0; t < ch.np; ++t) {
                    J[2 * m[t]] += f * Jc[2 * t];
                    J[2 * m[t] + 1] += f * Jc[2 * t + 1];
                }
            }
            m += ch.np;
        };
        if (n.kind == 2) {
            const double g = std::min(1.0 - 1e-9, std::max(1e-9, params[n.H]));
            addchild(n.child[0], g);
            addchild(n.child[1], 1 - g);
            const double *a = val + 2 * n.child[0], *b = val + 2 * n.child[1];
            J[2 * n.hpos] += a[0] - b[0];
            J[2 * n.hpos + 1] += a[1] - b[1];
        } else {
            const double early = std::max(n.set < 0 ? 1.0 : w.ex[static_cast<size_t>(n.set)], 1e-9);
            splitprobabilities(params[n.H], alpha, sp, sg);
            double below0 = 0.0, below1 = 0.0, dg0 = 0.0, dg1 = 0.0;
            for (int k = 0; k < 4; ++k) {
                const double f = early * sp[k];
                addchild(n.child[k], f);
                const double* c = val + 2 * n.child[k];
                below0 += f * c[0];
                below1 += f * c[1];
                dg0 += sg[k] * c[0];
                dg1 += sg[k] * c[1];
            }
            J[2 * n.hpos] += early * dg0;
            J[2 * n.hpos + 1] += early * dg1;
            if (n.set >= 0) {
                const uint32_t ns = setoff[static_cast<size_t>(n.set) + 1] - setoff[static_cast<size_t>(n.set)];
                const double h0 = (n.wc == 1 ? early : 0.0) - below0, h1 = (n.wc == 2 ? early : 0.0) - below1;
                for (uint32_t t = 0; t < ns; ++t) {
                    J[2 * m[t]] += h0;
                    J[2 * m[t] + 1] += h1;
                }
            }
        }
    }
}

}  // namespace csnaq
