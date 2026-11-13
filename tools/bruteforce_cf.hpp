// Independent expected quartet CFs of a network at ρ = 0, by brute force.
// Used to validate the quartet equations (tools/classcheck.cpp and the tests).
#pragma once

#include <array>
#include <cmath>
#include <cstdint>
#include <vector>

#include "network.hpp"

namespace csnaq {

// ------------------------------------------------------------------------
// Before the first coalescence the four lineages are distinct, and each picks
// a parent at every hybrid node independently (probability γ). For each
// combination of four root paths, edges are visited from the leaves up: in an
// edge shared by k >= 2 lineages over length t, a coalescence happens with
// probability 1 - exp(-C(k,2) t), between a pair chosen uniformly. The first
// coalescence decides the quartet topology; if none happens below the root,
// the three topologies are equally likely.
// ------------------------------------------------------------------------
struct BruteForceCF {
    const Network& net;
    std::vector<int> rank;            // edge -> position in bottom-up order
    std::vector<EdgeId> byrank;
    std::vector<double> npaths;       // node -> number of paths from the root
    explicit BruteForceCF(const Network& n) : net(n), rank(n.epool.size(), -1), npaths(n.npool.size(), -1.0) {
        // nodes in an order where parents come first
        std::vector<NodeId> order;
        std::vector<int> pending(net.npool.size(), 0);
        for (NodeId x : net.node)
            for (EdgeId e : net.N(x).edge)
                if (net.ischildof(x, e)) ++pending[static_cast<size_t>(x)];
        std::vector<NodeId> ready{net.getroot()};
        while (!ready.empty()) {
            NodeId x = ready.back();
            ready.pop_back();
            order.push_back(x);
            for (EdgeId e : net.N(x).edge)
                if (net.isparentof(x, e) && --pending[static_cast<size_t>(net.getchild(e))] == 0)
                    ready.push_back(net.getchild(e));
        }
        if (order.size() != net.node.size()) throw SnaqError("BruteForceCF: the network is not a rooted DAG");
        for (auto it = order.rbegin(); it != order.rend(); ++it)
            for (EdgeId e : net.N(*it).edge)
                if (net.ischildof(*it, e)) {
                    rank[static_cast<size_t>(e)] = static_cast<int>(byrank.size());
                    byrank.push_back(e);
                }
        for (NodeId x : order) {
            double c = x == net.getroot() ? 1.0 : 0.0;
            for (EdgeId e : net.N(x).edge)
                if (net.ischildof(x, e)) c += npaths[static_cast<size_t>(net.getparent(e))];
            npaths[static_cast<size_t>(x)] = c;
        }
    }
    struct Path {
        std::vector<uint64_t> bits;
        double prob;
    };
    void paths(NodeId x, std::vector<uint64_t>& cur, double prob, std::vector<Path>& out) const {
        bool any = false;
        for (EdgeId e : net.N(x).edge) {
            if (!net.ischildof(x, e)) continue;
            any = true;
            const size_t r = static_cast<size_t>(rank[static_cast<size_t>(e)]);
            cur[r >> 6] |= 1ULL << (r & 63);
            paths(net.getparent(e), cur, prob * (net.E(e).hybrid ? net.E(e).gamma : 1.0), out);
            cur[r >> 6] &= ~(1ULL << (r & 63));
        }
        if (!any) out.push_back({cur, prob});
    }
    // false if there are too many path combinations
    bool ecf(const std::array<NodeId, 4>& leaves, double* cf, double maxcombos = 3e5) const {
        double combos = 1.0;
        for (NodeId l : leaves) combos *= npaths[static_cast<size_t>(l)];
        if (combos > maxcombos) return false;
        const size_t W = (byrank.size() + 63) / 64;
        std::array<std::vector<Path>, 4> P;
        for (size_t k = 0; k < 4; ++k) {
            std::vector<uint64_t> cur(W, 0);
            paths(leaves[k], cur, 1.0, P[k]);
        }
        cf[0] = cf[1] = cf[2] = 0.0;
        static const int SPLIT[4][4] = {{-1, 0, 1, 2}, {0, -1, 2, 1}, {1, 2, -1, 0}, {2, 1, 0, -1}};
        for (const Path& a : P[0])
            for (const Path& b : P[1])
                for (const Path& c : P[2])
                    for (const Path& d : P[3]) {
                        const Path* ps[4] = {&a, &b, &c, &d};
                        double pnone = a.prob * b.prob * c.prob * d.prob;
                        for (size_t w = 0; w < W && pnone > 0; ++w) {
                            const uint64_t x0 = a.bits[w], x1 = b.bits[w], x2 = c.bits[w], x3 = d.bits[w];
                            uint64_t shared = (x0 & x1) | (x0 & x2) | (x0 & x3) | (x1 & x2) | (x1 & x3) | (x2 & x3);
                            while (shared) {
                                const size_t r = (w << 6) + static_cast<size_t>(ctz64(shared));
                                shared &= shared - 1;
                                int in[4], k = 0;
                                for (int j = 0; j < 4; ++j)
                                    if ((ps[j]->bits[r >> 6] >> (r & 63)) & 1) in[k++] = j;
                                const double npairs = k * (k - 1) / 2.0;
                                const double stay = std::exp(-npairs * std::max(net.E(byrank[r]).length, 0.0));
                                const double each = pnone * (1.0 - stay) / npairs;
                                for (int x = 0; x < k; ++x)
                                    for (int y = x + 1; y < k; ++y) cf[SPLIT[in[x]][in[y]]] += each;
                                pnone *= stay;
                            }
                        }
                        for (int k = 0; k < 3; ++k) cf[k] += pnone / 3.0;
                    }
        return true;
    }
};

}  // namespace csnaq
