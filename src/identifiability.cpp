#include "identifiability.hpp"

#include <algorithm>
#include <deque>
#include <string>
#include <tuple>

#include "network_manip.hpp"
#include "snaq_network.hpp"

namespace csnaq {

namespace {

bool contains(const std::vector<NodeId>& v, NodeId n) { return std::find(v.begin(), v.end(), n) != v.end(); }

// [bcp for bcp in biconnectedcomponents(N) if length(bcp) > 1]
std::vector<std::vector<EdgeId>> nontrivialblobs(Network& net) {
    std::vector<std::vector<EdgeId>> out;
    for (auto& comp : biconnectedcomponents(net, false))
        if (comp.size() > 1) out.push_back(std::move(comp));
    return out;
}

// attachednodes(n) = [e.node[1] == n ? e.node[2] : e.node[1] for e in n.edge]
std::vector<NodeId> attachednodes(const Network& net, NodeId n) {
    std::vector<NodeId> out;
    for (EdgeId e : net.N(n).edge) {
        const Edge& ed = net.E(e);
        out.push_back(ed.nd[0] == n ? ed.nd[1] : ed.nd[0]);
    }
    return out;
}

// PhyloNetworks.getlevel(net) (preorder = true, preprocess = true): directs the
// edges from the root (a side effect on `net`, as in Julia), then returns the
// maximum number of minor edges in a biconnected component.
int phylonetworks_getlevel(Network& net) {
    directedges(net);
    preorder(net);
    int level = 0;
    for (const auto& comp : biconnectedcomponents(net, false)) {
        int nminor = 0;
        for (EdgeId e : comp)
            if (!net.E(e).ismajor) ++nminor;
        level = std::max(level, nminor);
    }
    return level;
}

}  // namespace

std::vector<NodeId> getleavesbelow(const Network& net, const std::vector<NodeId>& start) {
    std::deque<NodeId> nodes(start.begin(), start.end());
    std::vector<NodeId> leaves;
    std::vector<NodeId> visited;
    while (!nodes.empty()) {
        NodeId curr = nodes.front();
        nodes.pop_front();
        if (contains(visited, curr)) continue;
        visited.push_back(curr);
        if (net.N(curr).leaf) {
            leaves.push_back(curr);
            continue;
        }
        for (NodeId c : net.getchildren(curr)) nodes.push_back(c);
    }
    return leaves;
}

std::vector<int> getdiamondns(const Network& net, const std::vector<NodeId>& cycle) {
    if (cycle.size() != 4) throw SnaqError("cycle must be of length 4 (has length " + std::to_string(cycle.size()) + ")");
    NodeId n3root = NONE, n1root = NONE;
    for (NodeId n : cycle)
        if (net.N(n).hybrid) {
            n3root = n;
            break;
        }
    if (n3root == NONE) throw SnaqError("getdiamondns: no hybrid node in the cycle");
    for (NodeId n : cycle) {
        if (n == n3root) continue;
        bool touches = false;
        for (EdgeId e : net.N(n).edge) {
            const Edge& ed = net.E(e);
            if (ed.nd[0] == n3root || ed.nd[1] == n3root) touches = true;
        }
        if (!touches) {
            n1root = n;
            break;
        }
    }
    if (n1root == NONE) throw SnaqError("getdiamondns: no node opposite the hybrid node");
    std::vector<NodeId> others;
    for (NodeId n : cycle)
        if (n != n3root && n != n1root) others.push_back(n);
    if (others.size() != 2) throw SnaqError("getdiamondns: malformed 4-cycle");
    const NodeId roots[4] = {others[0], n1root, others[1], n3root};

    std::vector<int> ns(4, 0);
    for (int idx = 0; idx < 4; ++idx) {
        // Same traversal as Julia: nodes are marked visited when dequeued, so a
        // node can be queued (and counted) more than once.
        std::vector<NodeId> visited(roots, roots + 4);
        std::deque<NodeId> Q{roots[idx]};
        while (!Q.empty()) {
            NodeId curr = Q.front();
            Q.pop_front();
            visited.push_back(curr);
            if (net.N(curr).leaf) {
                ns[static_cast<size_t>(idx)] += 1;
            } else {
                for (NodeId a : attachednodes(net, curr)) {
                    if (contains(visited, a)) continue;
                    Q.push_back(a);
                }
            }
        }
    }
    return ns;
}

std::pair<bool, int> isCk(const Network& input) {
    Network N = input.deepcopynetwork();
    if (!isgallednetwork(N)) return {false, 0};
    if (!std::get<2>(istreechild(N))) return {false, 0};
    if (shrink2cycles(N)) return {false, 0};

    int k = N.numtaxa;
    for (const auto& blob : nontrivialblobs(N)) {
        // lowest nodes of the blob = its hybrid nodes (N is galled)
        std::vector<NodeId> lnodes;
        for (EdgeId e : blob) {
            NodeId c = N.getchild(e);
            if (N.N(c).hybrid && !contains(lnodes, c)) lnodes.push_back(c);
        }
        std::vector<NodeId> Y_l = getleavesbelow(N, lnodes);
        for (size_t l = 0; l < Y_l.size(); ++l) {
            Network bloblet = N.deepcopynetwork();
            DeleteOpts opts;  // simplify = true, nofuse = false, multgammas = false
            for (size_t i = 0; i < Y_l.size(); ++i) {
                if (i == l) continue;
                deleteleaf_name(bloblet, N.N(Y_l[i]).name, opts);
            }
            auto bcp = nontrivialblobs(bloblet);
            if (bcp.size() != 1) return {false, 0};
            k = std::min(k, static_cast<int>(bcp[0].size()));
        }
    }
    return {true, k};
}

bool tcgidentifiable(const Network& input) {
    if (input.numhybrids == 0) return true;
    if (input.isrooted) return false;

    Network N = input.deepcopynetwork();
    if (shrink2cycles(N) || shrink3cycles(N)) return false;

    if (phylonetworks_getlevel(N) == 1) {
        // bad diamonds: a 4-cycle is good iff n0 >= 2 || n2 >= 2 || (n1 >= 2 && n3 >= 2)
        for (const auto& blob : nontrivialblobs(N)) {
            std::vector<NodeId> cycle;
            for (EdgeId e : blob) {
                for (NodeId n : {N.getchild(e), N.getparent(e)})
                    if (!contains(cycle, n)) cycle.push_back(n);
            }
            if (cycle.size() != 4) continue;
            auto ns = getdiamondns(N, cycle);
            if (!(ns[0] >= 2 || ns[2] >= 2 || (ns[1] >= 2 && ns[3] >= 2))) return false;
        }
    } else {
        auto [boolCk, k] = isCk(N);
        if (!boolCk) return false;
        if (k < 4) return false;
        if (k > 4) return true;
        // C4: need at least 2 taxa below each hybrid (lowest node of each blob)
        for (NodeId ln : N.hybrid)
            if (getleavesbelow(N, {ln}).size() < 2) return false;
    }
    return true;
}

}  // namespace csnaq
