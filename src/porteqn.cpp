#include "porteqn.hpp"

#include <algorithm>

namespace csnaq {

namespace {

using EdgeChain = SmallVec<EdgeId, 24>;

// Split (1: 01|23, 2: 02|13, 3: 03|12) that pairs lineages x and y.
const int PAIRSPLIT[4][4] = {{0, 1, 2, 3}, {1, 0, 3, 2}, {2, 3, 0, 1}, {3, 2, 1, 0}};

inline size_t commonsuffix(const EdgeChain& a, const EdgeChain& b) {
    size_t s = 0;
    const size_t na = a.size(), nb = b.size();
    while (s < na && s < nb && a[na - 1 - s] == b[nb - 1 - s]) ++s;
    return s;
}

}  // namespace

size_t equationsize(const RecursiveCFEquation& eqn) {
    size_t n = 1;
    for (const auto& d : eqn.divisions) n += equationsize(d);
    return n;
}

PortEquationBuilder::PortEquationBuilder(const Network& net) {
    const size_t NP = net.npool.size(), EP = net.epool.size(), N = net.node.size();
    pmaj_.assign(NP, NONE);
    pmin_.assign(NP, NONE);
    npid_.assign(NP, -1);
    npos_.assign(NP, -1);
    eparent_.assign(EP, NONE);
    echild_.assign(EP, NONE);
    epid_.assign(EP, -1);
    for (EdgeId e : net.edge) {
        eparent_[static_cast<size_t>(e)] = net.getparent(e);
        echild_[static_cast<size_t>(e)] = net.getchild(e);
        epid_[static_cast<size_t>(e)] = net.E(e).pid;
    }
    for (size_t j = 0; j < N; ++j) {
        const NodeId x = net.node[j];
        npos_[static_cast<size_t>(x)] = static_cast<int>(j);
        npid_[static_cast<size_t>(x)] = net.N(x).pid;
        int np = 0;
        EdgeId pe[2] = {NONE, NONE};
        for (EdgeId e : net.N(x).edge)
            if (net.ischildof(x, e)) {
                if (np < 2) pe[np] = e;
                ++np;
            }
        if (np > 2) throw SnaqError("PortEquationBuilder: a node has more than 2 parents");
        if (np == 1) {
            pmaj_[static_cast<size_t>(x)] = pe[0];
        } else if (np == 2) {
            const bool m0 = net.E(pe[0]).ismajor, m1 = net.E(pe[1]).ismajor;
            if (m0 == m1) throw SnaqError("PortEquationBuilder: a hybrid node needs one major and one minor parent edge");
            pmaj_[static_cast<size_t>(x)] = m0 ? pe[0] : pe[1];
            pmin_[static_cast<size_t>(x)] = m0 ? pe[1] : pe[0];
        }
    }
    // every node with its ancestors, parents first
    W_ = (N + 63) / 64;
    anc_.assign(N * W_, 0);
    std::vector<char> done(NP, 0);
    std::vector<NodeId> stack;
    for (NodeId start : net.node) {
        if (done[static_cast<size_t>(start)]) continue;
        stack.assign(1, start);
        while (!stack.empty()) {
            const NodeId x = stack.back();
            const size_t sx = static_cast<size_t>(x);
            if (done[sx]) {
                stack.pop_back();
                continue;
            }
            bool ready = true;
            for (EdgeId e : {pmaj_[sx], pmin_[sx]})
                if (e != NONE && !done[static_cast<size_t>(eparent_[static_cast<size_t>(e)])]) {
                    stack.push_back(eparent_[static_cast<size_t>(e)]);
                    ready = false;
                }
            if (!ready) {
                if (stack.size() > 4 * N + 16) throw SnaqError("PortEquationBuilder: the network has a directed cycle");
                continue;
            }
            const size_t bx = static_cast<size_t>(npos_[sx]);
            uint64_t* dst = &anc_[bx * W_];
            for (EdgeId e : {pmaj_[sx], pmin_[sx]})
                if (e != NONE) {
                    const uint64_t* src =
                        &anc_[static_cast<size_t>(npos_[static_cast<size_t>(eparent_[static_cast<size_t>(e)])]) * W_];
                    for (size_t w = 0; w < W_; ++w) dst[w] |= src[w];
                }
            dst[bx >> 6] |= 1ULL << (bx & 63);
            done[sx] = 1;
            stack.pop_back();
        }
    }
    maxdepth_ = static_cast<int>(8 * N + 64);
}

RecursiveCFEquation PortEquationBuilder::build(const std::array<NodeId, 4>& ports, NodeId top) const {
    State st;
    for (size_t k = 0; k < 4; ++k) {
        if (!isancestor(top, ancestors(Pos{ports[k], NONE})))
            throw SnaqError("PortEquationBuilder: `top` is not above a port");
        st[k] = Pos{ports[k], NONE};
    }
    return recur(st, top, 0);
}

RecursiveCFEquation PortEquationBuilder::recur(const State& start, NodeId top, int depth) const {
    State st = start;
    for (;; ++depth) {
        if (depth > maxdepth_) throw SnaqError("PortEquationBuilder: the recursion does not end");
        // ---- each lineage climbs to a hybrid node it has not passed, or to `top` ----
        EdgeChain chain[4];
        NodeId end[4];
        for (int k = 0; k < 4; ++k) {
            NodeId x = st[static_cast<size_t>(k)].node;
            if (st[static_cast<size_t>(k)].via != NONE) {
                chain[k].push_back(st[static_cast<size_t>(k)].via);
                x = eparent_[static_cast<size_t>(st[static_cast<size_t>(k)].via)];
            }
            while (x != top && pmin_[static_cast<size_t>(x)] == NONE) {
                const EdgeId e = pmaj_[static_cast<size_t>(x)];
                if (e == NONE) throw SnaqError("PortEquationBuilder: a lineage passed `top`");
                chain[k].push_back(e);
                x = eparent_[static_cast<size_t>(e)];
            }
            end[k] = x;
        }
        // ---- lineages that meet along tree-like paths: those with the same end ----
        int members[4], nmembers = 0, outsider = -1;
        for (int k = 0; k < 4 && nmembers < 3; ++k) {
            nmembers = 0;
            for (int j = 0; j < 4; ++j)
                if (end[j] == end[k]) members[nmembers++] = j;
        }
        if (nmembers == 3)
            for (int j = 0; j < 4; ++j)
                if (end[j] != end[members[0]]) outsider = j;
        if (nmembers >= 3) {
            // the two lineages that meet first, and the edges they share before a third one joins
            size_t s[4][4];
            int x = members[0], y = members[1];
            size_t best = 0, sall = ~size_t{0};
            for (int a = 0; a < nmembers; ++a)
                for (int b = a + 1; b < nmembers; ++b) {
                    const size_t v = commonsuffix(chain[members[a]], chain[members[b]]);
                    s[members[a]][members[b]] = s[members[b]][members[a]] = v;
                    sall = std::min(sall, v);
                    if (v > best) {
                        best = v;
                        x = members[a];
                        y = members[b];
                    }
                }
            bool tree = true;
            if (nmembers == 3) {
                // the fourth lineage must not be able to join the three before they have all met
                const uint64_t* ao = ancestors(st[static_cast<size_t>(outsider)]);
                for (int a = 0; a < 3 && tree; ++a) {
                    const EdgeChain& c = chain[members[a]];
                    for (size_t i = 0; i + sall < c.size() && tree; ++i)
                        if (isancestor(echild_[static_cast<size_t>(c[i])], ao)) tree = false;
                }
            }
            if (tree) {
                RecursiveCFEquation eq;
                eq.can_coalesce_here = true;
                eq.division_H = -1;
                eq.which_coal = PAIRSPLIT[x][y];
                // edges climbed by exactly two lineages: x and y, then the other two (if they are a pair)
                size_t lo = 0;
                for (int a = 0; a < nmembers; ++a)
                    if (members[a] != x && members[a] != y) lo = std::max(lo, s[x][members[a]]);
                const EdgeChain& cx = chain[x];
                for (size_t i = cx.size() - best; i + lo < cx.size(); ++i)
                    eq.coal_edges.push_back(epid_[static_cast<size_t>(cx[i])]);
                if (nmembers == 4) {
                    int z = -1, w = -1;
                    for (int a = 0; a < 4; ++a)
                        if (a != x && a != y) (z < 0 ? z : w) = a;
                    const size_t szw = s[z][w], lz = std::max(s[z][x], s[z][y]);
                    const EdgeChain& cz = chain[z];
                    for (size_t i = cz.size() - szw; i + lz < cz.size(); ++i)
                        eq.coal_edges.push_back(epid_[static_cast<size_t>(cz[i])]);
                }
                return eq;
            }
        }
        // ---- a lowest end hybrid node ----
        // (one with a single lineage if there is one: dividing in 2 before
        // dividing in 4 gives a smaller equation than the reverse)
        int first = -1, firstcount = 0;
        for (int k = 0; k < 4 && firstcount != 1; ++k) {
            if (end[k] == top) continue;
            bool lowest = true;
            int count = 0;
            for (int j = 0; j < 4 && lowest; ++j) {
                if (end[j] == end[k])
                    ++count;
                else if (end[j] != top &&
                         isancestor(end[k], &anc_[static_cast<size_t>(npos_[static_cast<size_t>(end[j])]) * W_]))
                    lowest = false;
            }
            if (lowest && (first < 0 || count < firstcount)) {
                first = k;
                firstcount = count;
            }
        }
        if (first < 0) throw SnaqError("PortEquationBuilder: no hybrid node to pass");
        const NodeId H = end[first];
        const EdgeId emin = pmin_[static_cast<size_t>(H)], emaj = pmaj_[static_cast<size_t>(H)];
        const NodeId umin = eparent_[static_cast<size_t>(emin)], umaj = eparent_[static_cast<size_t>(emaj)];
        int at[4], nat = 0;
        for (int k = 0; k < 4; ++k)
            if (end[k] == H) at[nat++] = k;
        if (nat == 1) {
            const int x = at[0];
            // The nodes where the lineage can first meet another one: if there
            // is only one, it gets there alone whatever the edges it takes.
            static thread_local std::vector<uint64_t> others, seen;
            static thread_local std::vector<NodeId> todo;
            others.assign(W_, 0);
            seen.assign(W_, 0);
            for (int k = 0; k < 4; ++k) {
                if (k == x) continue;
                const uint64_t* a = ancestors(st[static_cast<size_t>(k)]);
                for (size_t w = 0; w < W_; ++w) others[w] |= a[w];
            }
            todo.clear();
            todo.push_back(umin);
            todo.push_back(umaj);
            NodeId contact = NONE;
            bool single = true;
            while (!todo.empty() && single) {
                const NodeId r = todo.back();
                todo.pop_back();
                const size_t b = static_cast<size_t>(npos_[static_cast<size_t>(r)]);
                if ((seen[b >> 6] >> (b & 63)) & 1) continue;
                seen[b >> 6] |= 1ULL << (b & 63);
                if (isancestor(r, others.data())) {
                    if (contact == NONE)
                        contact = r;
                    else if (contact != r)
                        single = false;
                    continue;
                }
                for (EdgeId e : {pmaj_[static_cast<size_t>(r)], pmin_[static_cast<size_t>(r)]})
                    if (e != NONE) todo.push_back(eparent_[static_cast<size_t>(e)]);
            }
            if (contact == NONE) throw SnaqError("PortEquationBuilder: a lineage meets no other below `top`");
            if (single) {
                st[static_cast<size_t>(x)] = Pos{contact, NONE};
                continue;
            }
            RecursiveCFEquation eq;
            eq.can_coalesce_here = false;
            eq.which_coal = 0;
            eq.division_H = npid_[static_cast<size_t>(H)];
            eq.divisions.reserve(2);
            State div = st;
            div[static_cast<size_t>(x)] = Pos{umin, NONE};
            eq.divisions.push_back(recur(div, top, depth + 1));
            div[static_cast<size_t>(x)] = Pos{umaj, NONE};
            eq.divisions.push_back(recur(div, top, depth + 1));
            return eq;
        }
        if (nat != 2) throw SnaqError("PortEquationBuilder: unexpected number of lineages at a hybrid node");
        const int x = at[0], y = at[1];  // x < y
        RecursiveCFEquation eq;
        const size_t shared = commonsuffix(chain[x], chain[y]);
        for (size_t i = chain[x].size() - shared; i < chain[x].size(); ++i)
            eq.coal_edges.push_back(epid_[static_cast<size_t>(chain[x][i])]);
        eq.can_coalesce_here = shared > 0;
        eq.which_coal = PAIRSPLIT[x][y];
        eq.division_H = npid_[static_cast<size_t>(H)];
        eq.divisions.reserve(4);
        State div = st;
        div[static_cast<size_t>(x)] = div[static_cast<size_t>(y)] = Pos{H, emin};  // both take the minor edge
        eq.divisions.push_back(recur(div, top, depth + 1));
        div[static_cast<size_t>(x)] = div[static_cast<size_t>(y)] = Pos{H, emaj};  // both take the major edge
        eq.divisions.push_back(recur(div, top, depth + 1));
        div[static_cast<size_t>(x)] = Pos{umin, NONE};  // one each
        div[static_cast<size_t>(y)] = Pos{umaj, NONE};
        eq.divisions.push_back(recur(div, top, depth + 1));
        div[static_cast<size_t>(x)] = Pos{umaj, NONE};
        div[static_cast<size_t>(y)] = Pos{umin, NONE};
        eq.divisions.push_back(recur(div, top, depth + 1));
        return eq;
    }
}

}  // namespace csnaq
