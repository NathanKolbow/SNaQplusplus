#include "classes.hpp"

#include <algorithm>
#include <cmath>

#include "parallel.hpp"
#include "quartets.hpp"

namespace csnaq {

QuartetList QuartetList::allquartets(int ntaxa) {
    QuartetList q;
    q.ntaxa = ntaxa;
    q.all = true;
    return q;
}

QuartetList QuartetList::fromranks(int ntaxa, const std::vector<long long>& ranks) {
    QuartetList q;
    q.ntaxa = ntaxa;
    if (static_cast<long long>(ranks.size()) == nchoose4(ntaxa)) {
        bool iota = true;
        for (size_t j = 0; j < ranks.size() && iota; ++j) iota = ranks[j] == static_cast<long long>(j);
        if (iota) {
            q.all = true;
            return q;
        }
    }
    q.taxa.resize(ranks.size());
    for (size_t j = 0; j < ranks.size(); ++j) q.taxa[j] = unrank4taxa(ntaxa, ranks[j]);
    return q;
}

double obsqlogq(const CFMatrix& obs, const double* weights) {
    double total = 0.0;
    for (size_t j = 0; j < obs.nrow(); ++j) {
        double s = 0.0;
        for (int k = 0; k < 3; ++k) {
            const double q = obs(j, k);
            if (q > 0) s += q * std::log(q);
        }
        total += weights ? weights[j] * s : s;
    }
    return total;
}

namespace {

inline uint64_t mix64(uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}
inline uint64_t hashrec(uint64_t seed, uint64_t a, uint64_t b, uint64_t c, uint64_t d) {
    uint64_t h = mix64(seed ^ a);
    h = mix64(h ^ (b + 0x632be59bd9b4e019ULL));
    h = mix64(h ^ (c + 0x2545f4914f6cdd1dULL));
    return mix64(h ^ (d + 0x9fb21c651e98df25ULL));
}

void collectparams(const RecursiveCFEquation& e, std::vector<int>& out) {
    for (int p : e.coal_edges)
        if (std::find(out.begin(), out.end(), p) == out.end()) out.push_back(p);
    if (e.division_H >= 0 && std::find(out.begin(), out.end(), e.division_H) == out.end()) out.push_back(e.division_H);
    for (const auto& d : e.divisions) collectparams(d, out);
}

// Observed-CF column (0: 12|34, 1: 13|24, 2: 14|23) that pairs taxa x and y.
const uint8_t PAIRCOL[4][4] = {{0, 0, 1, 2}, {0, 0, 2, 1}, {1, 2, 0, 0}, {2, 1, 0, 0}};

// Open-addressing table: packed ports -> class index.
struct PortTable {
    std::vector<uint64_t> keys;
    std::vector<uint32_t> vals;
    size_t mask = 0, count = 0;
    static constexpr uint64_t EMPTY = ~0ULL;
    explicit PortTable(size_t cap = 1024) { reset(cap); }
    void reset(size_t cap) {
        keys.assign(cap, EMPTY);
        vals.assign(cap, 0);
        mask = cap - 1;
        count = 0;
    }
    // Index of `key`, inserted with value `next` if absent.
    uint32_t get(uint64_t key, uint32_t next, bool& isnew) {
        size_t i = static_cast<size_t>(mix64(key)) & mask;
        while (keys[i] != EMPTY) {
            if (keys[i] == key) {
                isnew = false;
                return vals[i];
            }
            i = (i + 1) & mask;
        }
        keys[i] = key;
        vals[i] = next;
        isnew = true;
        if (++count * 2 > keys.size()) grow();
        return next;
    }
    void grow() {
        std::vector<uint64_t> ok = std::move(keys);
        std::vector<uint32_t> ov = std::move(vals);
        reset(ok.size() * 2);
        for (size_t j = 0; j < ok.size(); ++j) {
            if (ok[j] == EMPTY) continue;
            size_t i = static_cast<size_t>(mix64(ok[j])) & mask;
            while (keys[i] != EMPTY) i = (i + 1) & mask;
            keys[i] = ok[j];
            vals[i] = ov[j];
            ++count;
        }
    }
};

}  // namespace

QuartetClasses buildclasses(Network& net, const QuartetList& ql, const CFMatrix& obs, const ClassOptions& opt) {
    gatheroptimizationinfo(net, true);  // unique numbers
    const size_t NP = net.npool.size(), EP = net.epool.size();
    const int n = ql.ntaxa;
    const size_t sn = static_cast<size_t>(n);
    if (obs.nrow() != ql.size()) throw SnaqError("buildclasses: one row of observed CFs is needed per quartet");

    std::vector<int> nodepos(NP, -1), edgepos(EP, -1);
    for (size_t j = 0; j < net.node.size(); ++j) nodepos[static_cast<size_t>(net.node[j])] = static_cast<int>(j);
    for (size_t j = 0; j < net.edge.size(); ++j) edgepos[static_cast<size_t>(net.edge[j])] = static_cast<int>(j);

    // taxa, sorted by name
    std::vector<std::pair<std::string, NodeId>> lv;
    for (NodeId l : net.leaf) lv.emplace_back(net.N(l).name, l);
    std::sort(lv.begin(), lv.end());
    if (static_cast<int>(lv.size()) != n) throw SnaqError("buildclasses: the network and the quartets differ in taxa");

    // ---- cut edges (bridges) of the underlying undirected graph ----
    std::vector<char> isbridge(EP, 0);
    {
        std::vector<int> disc(NP, -1), low(NP, 0);
        struct Frame {
            NodeId node;
            EdgeId pe;
            size_t next;
        };
        std::vector<Frame> st;
        int t = 0;
        NodeId root = net.getroot();
        disc[static_cast<size_t>(root)] = low[static_cast<size_t>(root)] = t++;
        st.push_back({root, NONE, 0});
        while (!st.empty()) {
            Frame& f = st.back();
            const auto& edges = net.N(f.node).edge;
            if (f.next < edges.size()) {
                EdgeId e = edges[f.next++];
                if (e == f.pe) continue;
                NodeId o = net.getOtherNode(e, f.node);
                if (disc[static_cast<size_t>(o)] < 0) {
                    disc[static_cast<size_t>(o)] = low[static_cast<size_t>(o)] = t++;
                    st.push_back({o, e, 0});
                } else {
                    low[static_cast<size_t>(f.node)] =
                        std::min(low[static_cast<size_t>(f.node)], disc[static_cast<size_t>(o)]);
                }
            } else {
                Frame done = f;
                st.pop_back();
                if (!st.empty()) {
                    NodeId par = st.back().node;
                    low[static_cast<size_t>(par)] =
                        std::min(low[static_cast<size_t>(par)], low[static_cast<size_t>(done.node)]);
                    if (low[static_cast<size_t>(done.node)] > disc[static_cast<size_t>(par)])
                        isbridge[static_cast<size_t>(done.pe)] = 1;
                }
            }
        }
    }
    std::vector<EdgeId> bridges;
    for (EdgeId e : net.edge)
        if (isbridge[static_cast<size_t>(e)]) bridges.push_back(e);
    std::sort(bridges.begin(), bridges.end(), [&](EdgeId x, EdgeId y) { return net.E(x).pid < net.E(y).pid; });
    if (bridges.size() >= 65535) throw SnaqError("buildclasses: too many edges");
    std::vector<uint16_t> brank(EP, 0);
    for (size_t r = 0; r < bridges.size(); ++r) brank[static_cast<size_t>(bridges[r])] = static_cast<uint16_t>(r);

    // ---- the cut edges above each taxon, from the root down ----
    std::vector<EdgeId> anype(NP, NONE);
    for (EdgeId e : net.edge) anype[static_cast<size_t>(net.getchild(e))] = e;
    std::vector<uint32_t> boff(sn + 1, 0);
    std::vector<uint16_t> blist;
    for (int t = 0; t < n; ++t) {
        size_t start = blist.size();
        NodeId node = lv[static_cast<size_t>(t)].second;
        size_t guard = 0;
        while (anype[static_cast<size_t>(node)] != NONE) {
            EdgeId e = anype[static_cast<size_t>(node)];
            if (isbridge[static_cast<size_t>(e)]) blist.push_back(brank[static_cast<size_t>(e)]);
            node = net.getparent(e);
            if (++guard > NP) throw SnaqError("buildclasses: the network has a directed cycle");
        }
        std::reverse(blist.begin() + static_cast<std::ptrdiff_t>(start), blist.end());
        boff[static_cast<size_t>(t) + 1] = static_cast<uint32_t>(blist.size());
    }
    // lcp[a][b]: number of cut edges above both a and b
    std::vector<uint16_t> lcp(sn * sn, 0);
    for (size_t a = 0; a < sn; ++a)
        for (size_t b = a + 1; b < sn; ++b) {
            const uint16_t* pa = &blist[boff[a]];
            const uint16_t* pb = &blist[boff[b]];
            const size_t m = std::min(boff[a + 1] - boff[a], boff[b + 1] - boff[b]);
            size_t k = 0;
            while (k < m && pa[k] == pb[k]) ++k;
            lcp[a * sn + b] = lcp[b * sn + a] = static_cast<uint16_t>(k);
        }

    // ---- one pass over the quartets ----
    QuartetClasses out;
    out.npid = net.next_pid;
    PortTable table;
    std::vector<uint64_t> classports;
    std::vector<uint16_t> classtop;  // lowest cut edge above the four ports (NOTOP: none)
    const uint16_t NOTOP = 0xffff;
    const double* weights = opt.weights;
    if (opt.perquartet) {
        out.qclass.resize(ql.size());
        out.qcols.resize(ql.size());
    }
    const uint16_t* bl = blist.data();
    const uint32_t* bo = boff.data();
    const uint16_t* lc = lcp.data();
    auto addquartet = [&](size_t j, size_t a, size_t b, size_t c, size_t d) {
        const uint16_t lab = lc[a * sn + b], lac = lc[a * sn + c], lad = lc[a * sn + d], lbc = lc[b * sn + c],
                       lbd = lc[b * sn + d], lcd = lc[c * sn + d];
        uint16_t pr[4] = {bl[bo[a] + std::max({lab, lac, lad})], bl[bo[b] + std::max({lab, lbc, lbd})],
                          bl[bo[c] + std::max({lac, lbc, lcd})], bl[bo[d] + std::max({lad, lbd, lcd})]};
        int o[4] = {0, 1, 2, 3};
        for (int x = 1; x < 4; ++x)
            for (int y = x; y > 0 && pr[o[y - 1]] > pr[o[y]]; --y) std::swap(o[y - 1], o[y]);
        const uint64_t key = static_cast<uint64_t>(pr[o[0]]) | (static_cast<uint64_t>(pr[o[1]]) << 16) |
                             (static_cast<uint64_t>(pr[o[2]]) << 32) | (static_cast<uint64_t>(pr[o[3]]) << 48);
        bool isnew = false;
        const uint32_t idx = table.get(key, static_cast<uint32_t>(classports.size()), isnew);
        if (isnew) {
            classports.push_back(key);
            out.Q.push_back({0.0, 0.0, 0.0});
            const uint16_t m = std::min({lab, lac, lad, lbc, lbd, lcd});
            classtop.push_back(m > 0 ? bl[bo[a] + m - 1] : NOTOP);
        }
        const uint8_t c0 = PAIRCOL[o[0]][o[1]], c1 = PAIRCOL[o[0]][o[2]], c2 = PAIRCOL[o[0]][o[3]];
        const double* q = obs.row(j);
        auto& Q = out.Q[idx];
        if (weights) {
            const double w = weights[j];
            if (w != 0.0) {
                Q[0] += w * q[c0];
                Q[1] += w * q[c1];
                Q[2] += w * q[c2];
            }
        } else {
            Q[0] += q[c0];
            Q[1] += q[c1];
            Q[2] += q[c2];
        }
        if (opt.perquartet) {
            out.qclass[j] = idx;
            out.qcols[j] = static_cast<uint8_t>(c0 | (c1 << 2) | (c2 << 4));
        }
    };
    if (ql.all) {
        size_t j = 0;
        for (size_t d = 3; d < sn; ++d)
            for (size_t c = 2; c < d; ++c)
                for (size_t b = 1; b < c; ++b)
                    for (size_t a = 0; a < b; ++a) addquartet(j++, a, b, c, d);
    } else {
        for (size_t j = 0; j < ql.taxa.size(); ++j) {
            const auto& t = ql.taxa[j];
            addquartet(j, static_cast<size_t>(t[0]), static_cast<size_t>(t[1]), static_cast<size_t>(t[2]),
                       static_cast<size_t>(t[3]));
        }
    }
    const size_t nc = classports.size();

    // ---- the key of each class: its ports and every edge above them ----
    const size_t E = net.edge.size();
    const size_t W = (E + 63) / 64;
    std::vector<uint64_t> h1(E), h2(E);
    for (size_t pos = 0; pos < E; ++pos) {
        const Edge& ed = net.E(net.edge[pos]);
        const NodeId ch = net.getchild(net.edge[pos]), pa = net.getparent(net.edge[pos]);
        // (Edge::containroot is not part of the key: building an equation only
        // updates it, and never reads it to decide anything.)
        const uint64_t flags = (ed.hybrid ? 1u : 0u) | (ed.ismajor ? 2u : 0u) | (net.N(ch).hybrid ? 8u : 0u);
        const uint64_t pe = static_cast<uint64_t>(ed.pid), pp = static_cast<uint64_t>(net.N(pa).pid),
                       pc = static_cast<uint64_t>(net.N(ch).pid);
        h1[pos] = hashrec(0x1234567ULL, pe, pp, pc, flags);
        h2[pos] = hashrec(0xabcdef01ULL, pe, pp, pc, flags);
    }
    std::vector<int> ancidx(NP, -1);
    std::vector<uint64_t> ancstore;
    std::vector<NodeId> order;  // nodes whose ancestors are being computed (explicit stack)
    auto ancestors = [&](NodeId start) {
        if (ancidx[static_cast<size_t>(start)] >= 0) return;
        order.assign(1, start);
        while (!order.empty()) {
            NodeId x = order.back();
            if (ancidx[static_cast<size_t>(x)] >= 0) {
                order.pop_back();
                continue;
            }
            bool ready = true;
            for (EdgeId e : net.N(x).edge)
                if (net.ischildof(x, e) && ancidx[static_cast<size_t>(net.getparent(e))] < 0) {
                    order.push_back(net.getparent(e));
                    ready = false;
                }
            if (!ready) {
                if (order.size() > 4 * NP + 16) throw SnaqError("buildclasses: the network has a directed cycle");
                continue;
            }
            const size_t idx = ancstore.size() / std::max<size_t>(W, 1);
            ancstore.resize(ancstore.size() + W, 0);
            uint64_t* dst = &ancstore[idx * W];
            for (EdgeId e : net.N(x).edge)
                if (net.ischildof(x, e)) {
                    const uint64_t* src = &ancstore[static_cast<size_t>(ancidx[static_cast<size_t>(net.getparent(e))]) * W];
                    for (size_t w = 0; w < W; ++w) dst[w] |= src[w];
                    const size_t pos = static_cast<size_t>(edgepos[static_cast<size_t>(e)]);
                    dst[pos >> 6] |= 1ULL << (pos & 63);
                }
            ancidx[static_cast<size_t>(x)] = static_cast<int>(idx);
            order.pop_back();
        }
    };
    out.key.resize(nc);
    out.eqn.resize(nc);
    std::vector<std::array<NodeId, 4>> portnodes(nc);
    std::vector<NodeId> topnode(nc);
    std::vector<uint32_t> tobuild;
    std::vector<uint64_t> un(W);
    const NodeId rootnode = net.getroot();
    for (size_t ci = 0; ci < nc; ++ci) {
        ClassKey k;
        std::fill(un.begin(), un.end(), 0);
        // the equation depends on the network between `top` and the ports only
        const NodeId top = classtop[ci] == NOTOP ? rootnode : net.getchild(bridges[classtop[ci]]);
        topnode[ci] = top;
        ancestors(top);
        for (int x = 0; x < 4; ++x) {
            const EdgeId be = bridges[(classports[ci] >> (16 * x)) & 0xffff];
            const NodeId pn = net.getchild(be);
            portnodes[ci][static_cast<size_t>(x)] = pn;
            ancestors(pn);
            const uint64_t* src = &ancstore[static_cast<size_t>(ancidx[static_cast<size_t>(pn)]) * W];
            for (size_t w = 0; w < W; ++w) un[w] |= src[w];
            const uint64_t pp = static_cast<uint64_t>(net.E(be).pid);
            k.a += hashrec(0x51ed270bULL, static_cast<uint64_t>(x), pp, 0, 0);
            k.b ^= hashrec(0x2f9be6ccULL, static_cast<uint64_t>(x), pp, 0, 0);
        }
        const uint64_t* above = &ancstore[static_cast<size_t>(ancidx[static_cast<size_t>(top)]) * W];
        for (size_t w = 0; w < W; ++w) {
            uint64_t bits = un[w] & ~above[w];
            while (bits) {
                const size_t pos = (w << 6) + static_cast<size_t>(__builtin_ctzll(bits));
                bits &= bits - 1;
                k.a += h1[pos];
                k.b ^= h2[pos];
            }
        }
        out.key[ci] = k;
        if (opt.reuse) {
            auto it = opt.reuse->bykey.find(k);
            if (it != opt.reuse->bykey.end()) {
                out.eqn[ci] = opt.reuse->eqn[it->second];
                continue;
            }
        }
        tobuild.push_back(static_cast<uint32_t>(ci));
    }

    // ---- equations of the classes not found in `reuse` ----
    out.nbuilt = tobuild.size();
    if (!tobuild.empty()) {
        auto pend = std::make_shared<QuartetClasses::Pending>();
        for (EdgeId e : net.edge) pend->pm.set(net.E(e).number, net.E(e).pid);
        for (NodeId h : net.hybrid) pend->pm.set(net.N(h).number, net.N(h).pid);
        std::vector<int> chainidx(NP, -1);
        std::vector<Chain> chains;
        for (uint32_t ci : tobuild)
            for (NodeId pn : portnodes[ci])
                if (chainidx[static_cast<size_t>(pn)] < 0) {
                    chainidx[static_cast<size_t>(pn)] = static_cast<int>(chains.size());
                    chains.push_back(treelikechain(net, pn));
                }
        // tree-like classes are cheap; the others are built by finishclasses
        for (uint32_t ci : tobuild) {
            std::array<const Chain*, 4> ch;
            for (size_t x = 0; x < 4; ++x)
                ch[x] = &chains[static_cast<size_t>(chainidx[static_cast<size_t>(portnodes[ci][x])])];
            auto ce = std::make_shared<ClassEqn>();
            if (porttreelikeequation(net, ch, pend->pm, ce->eqn)) {
                collectparams(ce->eqn, ce->params);
                out.eqn[ci] = std::move(ce);
            } else {
                pend->classes.push_back(ci);
                pend->ports.push_back(portnodes[ci]);
                pend->top.push_back(topnode[ci]);
            }
        }
        if (!pend->classes.empty()) {
            pend->nodepos = std::move(nodepos);
            pend->edgepos = std::move(edgepos);
            out.pending = std::move(pend);
        }
    }
    out.batch.npid = out.npid;
    out.batch.qlogq = opt.qlogq;
    if (!opt.defer) finishclasses(net, out, opt.nthreads);
    return out;
}

void finishclasses(const Network& net, QuartetClasses& cls, int nthreads) {
    if (cls.pending) {
        const QuartetClasses::Pending& pend = *cls.pending;
        parallel_for(pend.classes.size(), nthreads, [&](size_t lo, size_t hi, int) {
            for (size_t t = lo; t < hi; ++t) {
                auto ce = std::make_shared<ClassEqn>();
                ce->eqn = portreticulateequation(net, pend.nodepos, pend.edgepos, pend.ports[t], pend.pm, pend.top[t]);
                collectparams(ce->eqn, ce->params);
                cls.eqn[pend.classes[t]] = std::move(ce);
            }
        });
        cls.pending.reset();
    }
    const size_t nc = cls.size();
    cls.bykey.clear();
    cls.bykey.reserve(nc * 2);
    for (size_t ci = 0; ci < nc; ++ci) cls.bykey.emplace(cls.key[ci], static_cast<uint32_t>(ci));

    // ---- batch for the fit ----
    ClassBatch& b = cls.batch;
    b.coal.clear();
    b.off.assign(1, 0);
    b.qres.clear();
    b.qoth.clear();
    b.gen.clear();
    b.genQ.clear();
    for (size_t ci = 0; ci < nc; ++ci) {
        const RecursiveCFEquation& eq = cls.eqn[ci]->eqn;
        const auto& Q = cls.Q[ci];
        if (eq.division_H == -1) {
            for (int p : eq.coal_edges) b.coal.push_back(p);
            b.off.push_back(static_cast<uint32_t>(b.coal.size()));
            const double r = Q[static_cast<size_t>(eq.which_coal - 1)];
            b.qres.push_back(r);
            b.qoth.push_back((Q[0] + Q[1] + Q[2]) - r);
        } else {
            b.gen.push_back(cls.eqn[ci].get());
            b.genQ.push_back(Q);
        }
    }
}

double classupperbound(const QuartetClasses& cls) {
    double ub = -cls.batch.qlogq;
    for (size_t ci = 0; ci < cls.size(); ++ci) {
        const auto& Q = cls.Q[ci];
        const ClassEqn* ce = cls.eqn[ci].get();
        if (ce && ce->eqn.division_H == -1) {
            const double r = Q[static_cast<size_t>(ce->eqn.which_coal - 1)];
            ub += treeclassbound(r, (Q[0] + Q[1] + Q[2]) - r, !ce->eqn.coal_edges.empty());
        } else {
            ub += generalclassbound(Q);
        }
    }
    return ub;
}

std::vector<double> pidparams(const Network& net, int npid) {
    std::vector<double> v(static_cast<size_t>(npid), 0.0);
    for (EdgeId e : net.edge) v[static_cast<size_t>(net.E(e).pid)] = net.E(e).length;
    for (NodeId h : net.hybrid) v[static_cast<size_t>(net.N(h).pid)] = net.E(net.getparentedgeminor(h)).gamma;
    return v;
}

double scoreclasses(const Network& net, const QuartetClasses& cls, double rho, int nthreads) {
    std::vector<double> params = pidparams(net, cls.npid), grad(static_cast<size_t>(cls.npid));
    LossWorkspace ws;
    return evalclasses(cls.batch, params.data(), grad.data(), rhotoalpha(rho), ws, nthreads);
}

}  // namespace csnaq
