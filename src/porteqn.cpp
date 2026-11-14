#include "porteqn.hpp"

#include <algorithm>

namespace csnaq {

namespace {

// Split (1: 01|23, 2: 02|13, 3: 03|12) that pairs lineages x and y.
const int PAIRSPLIT[4][4] = {{0, 1, 2, 3}, {1, 0, 3, 2}, {2, 3, 0, 1}, {3, 2, 1, 0}};

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
    tend_.assign(NP, NONE);
    bypos_.assign(N, NONE);
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
        bypos_[j] = x;
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
            // (its parents are done)
            tend_[sx] = (pmin_[sx] != NONE || pmaj_[sx] == NONE)
                            ? x
                            : tend_[static_cast<size_t>(eparent_[static_cast<size_t>(pmaj_[sx])])];
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
    const size_t W = W_;
    // P + k * W: the nodes of the tree-like path of lineage k, its end included
    SmallVec<uint64_t, 32> pathbuf;
    for (size_t i = 0; i < 4 * W; ++i) pathbuf.push_back(0);
    uint64_t* const P = pathbuf.data();
    const NodeId topend = tend_[static_cast<size_t>(top)];
    for (;; ++depth) {
        if (depth > maxdepth_) throw SnaqError("PortEquationBuilder: the recursion does not end");
        // ---- each lineage climbs to a hybrid node it has not passed, or to `top` ----
        NodeId end[4];
        size_t len[4];  // number of edges it climbs
        for (size_t k = 0; k < 4; ++k) {
            const Pos& p = st[k];
            const NodeId t = p.via != NONE ? eparent_[static_cast<size_t>(p.via)] : p.node;
            NodeId e = tend_[static_cast<size_t>(t)];
            if (e == topend)
                e = top;  // `top` is on the way
            else if (pmin_[static_cast<size_t>(e)] == NONE)
                throw SnaqError("PortEquationBuilder: a lineage passed `top`");
            end[k] = e;
            // the ancestors of t that are not above e
            const uint64_t* at = &anc_[static_cast<size_t>(npos_[static_cast<size_t>(t)]) * W];
            const size_t be = static_cast<size_t>(npos_[static_cast<size_t>(e)]);
            const uint64_t* ae = &anc_[be * W];
            uint64_t* dst = P + k * W;
            size_t n = 0;
            for (size_t w = 0; w < W; ++w) {
                dst[w] = at[w] & ~ae[w];
                n += static_cast<size_t>(popcount64(dst[w]));
            }
            dst[be >> 6] |= 1ULL << (be & 63);
            len[k] = n + (p.via != NONE ? 1 : 0);
        }
        // whether lineages x and y climb the same hybrid edge first (their paths are then the same)
        auto samevia = [&](int x, int y) {
            return st[static_cast<size_t>(x)].via != NONE &&
                   st[static_cast<size_t>(x)].via == st[static_cast<size_t>(y)].via;
        };
        // number of edges that lineages x and y, with the same end, both climb
        auto shared = [&](int x, int y) {
            const uint64_t *px = P + static_cast<size_t>(x) * W, *py = P + static_cast<size_t>(y) * W;
            size_t n = 0;
            for (size_t w = 0; w < W; ++w) n += static_cast<size_t>(popcount64(px[w] & py[w]));
            return n - 1 + (samevia(x, y) ? 1 : 0);
        };
        // the edges that x and y both climb, except those that `u` or `v` (if >= 0) climbs too
        auto sharededges = [&](int x, int y, int u, int v, bool withvia, ParamList& out) {
            const uint64_t *px = P + static_cast<size_t>(x) * W, *py = P + static_cast<size_t>(y) * W;
            const size_t be = static_cast<size_t>(npos_[static_cast<size_t>(end[x])]);
            if (withvia) out.push_back(epid_[static_cast<size_t>(st[static_cast<size_t>(x)].via)]);
            for (size_t w = 0; w < W; ++w) {
                uint64_t bits = px[w] & py[w];
                if (u >= 0) bits &= ~P[static_cast<size_t>(u) * W + w];
                if (v >= 0) bits &= ~P[static_cast<size_t>(v) * W + w];
                if (w == (be >> 6)) bits &= ~(1ULL << (be & 63));  // the end is not the child of an edge climbed
                while (bits) {
                    const size_t pos = (w << 6) + static_cast<size_t>(ctz64(bits));
                    bits &= bits - 1;
                    out.push_back(epid_[static_cast<size_t>(pmaj_[static_cast<size_t>(bypos_[pos])])]);
                }
            }
        };
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
                    const size_t v = shared(members[a], members[b]);
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
                // the fourth lineage must not be able to join the three before they have all met:
                // no node below the point where they have all met is one of its ancestors
                const uint64_t* ao = ancestors(st[static_cast<size_t>(outsider)]);
                const uint64_t *p0 = P + static_cast<size_t>(members[0]) * W,
                               *p1 = P + static_cast<size_t>(members[1]) * W,
                               *p2 = P + static_cast<size_t>(members[2]) * W;
                for (int a = 0; a < 3 && tree; ++a) {
                    const size_t m = static_cast<size_t>(members[a]);
                    const uint64_t* pa = P + m * W;
                    for (size_t w = 0; w < W && tree; ++w)
                        if (pa[w] & ~(p0[w] & p1[w] & p2[w]) & ao[w]) tree = false;
                    // (the hybrid node below the hybrid edge it climbs first)
                    if (tree && st[m].via != NONE && sall < len[m] && isancestor(st[m].node, ao)) tree = false;
                }
            }
            if (tree) {
                RecursiveCFEquation eq;
                eq.can_coalesce_here = true;
                eq.division_H = -1;
                eq.which_coal = PAIRSPLIT[x][y];
                // edges climbed by exactly two lineages: x and y, then the other two (if they are a pair)
                int z = -1, w = -1;
                for (int a = 0; a < nmembers; ++a)
                    if (members[a] != x && members[a] != y) (z < 0 ? z : w) = members[a];
                size_t lo = s[x][z];
                if (w >= 0) lo = std::max(lo, s[x][w]);
                sharededges(x, y, z, w, samevia(x, y) && lo < best, eq.coal_edges);
                if (nmembers == 4) {
                    const size_t szw = s[z][w], lz = std::max(s[z][x], s[z][y]);
                    sharededges(z, w, x, y, samevia(z, w) && lz < szw, eq.coal_edges);
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
            SmallVec<uint64_t, 8> others, seen;
            SmallVec<NodeId, 32> todo;
            for (size_t w = 0; w < W_; ++w) {
                others.push_back(0);
                seen.push_back(0);
            }
            for (int k = 0; k < 4; ++k) {
                if (k == x) continue;
                const uint64_t* a = ancestors(st[static_cast<size_t>(k)]);
                for (size_t w = 0; w < W_; ++w) others[w] |= a[w];
            }
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
        sharededges(x, y, -1, -1, samevia(x, y), eq.coal_edges);  // the edges they share below H
        eq.can_coalesce_here = !eq.coal_edges.empty();
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
