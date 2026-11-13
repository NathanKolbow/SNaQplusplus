#include "snaq_network.hpp"

#include <algorithm>
#include <map>
#include <set>

#include "network_manip.hpp"

namespace csnaq {

void semidirectnetwork(Network& net, bool check_conditions) {
    if (net.N(net.getroot()).edge.size() == 2) fuseedgesat(net, net.rooti, false);
    net.isrooted = false;
    for (EdgeId e : net.edge)
        if (!net.N(net.getchild(e)).leaf) net.E(e).length = std::max(net.E(e).length, 0.0);
    if (check_conditions) {
        for (NodeId n : net.node) {
            size_t deg = net.N(n).edge.size();
            if (net.N(n).leaf) {
                if (deg != 1) throw SnaqError("Found leaf in `N` with degree " + std::to_string(deg));
            } else if (deg != 3) {
                throw SnaqError("Found internal node in `N` with degrees " + std::to_string(deg));
            }
        }
        for (EdgeId e : net.edge)
            if (net.E(e).nnode != 2) throw SnaqError("Found edge with ≠2 attached nodes.");
    }
}

bool isdescendantof(const Network& net, NodeId descendant, NodeId ancestor) {
    std::vector<NodeId> queue{ancestor};
    long iters = 0;
    while (!queue.empty()) {
        NodeId curr = queue.back();
        queue.pop_back();
        if (curr == descendant) return true;
        for (NodeId c : net.getchildren(curr)) queue.push_back(c);
        if (++iters >= 100000) throw SnaqError("Looped " + std::to_string(iters) + " times in `isdescendantof`");
    }
    return false;
}

void swapedgehybridinfo(Network& net, EdgeId e1, EdgeId e2) {
    Edge& a = net.E(e1);
    Edge& b = net.E(e2);
    std::swap(a.gamma, b.gamma);
    std::swap(a.hybrid, b.hybrid);
    std::swap(a.ismajor, b.ismajor);
}

bool areconnected(const Network& net, NodeId n1, NodeId n2) {
    for (EdgeId e : net.N(n1).edge) {
        const Edge& ed = net.E(e);
        if (ed.nnode == 2 && ((ed.nd[0] == n1 && ed.nd[1] == n2) || (ed.nd[0] == n2 && ed.nd[1] == n1))) return true;
    }
    return false;
}

int getnetworklevel(Network& net) {
    if (net.numhybrids == 0) return 0;
    auto bcc = biconnectedcomponents(net);
    int max_level = 1;
    for (const auto& comp : bcc) {
        if (comp.size() <= 1) continue;
        std::vector<NodeId> seen;
        for (EdgeId e : comp) {
            NodeId child = net.getchild(e);
            if (net.N(child).hybrid && std::find(seen.begin(), seen.end(), child) == seen.end())
                seen.push_back(child);
        }
        max_level = std::max(max_level, static_cast<int>(seen.size()));
    }
    return max_level;
}

bool isgallednetwork(Network& net) {
    auto bcc = biconnectedcomponents(net);
    std::map<NodeId, int> counts;
    for (const auto& comp : bcc) {
        if (comp.size() <= 1) continue;
        for (EdgeId e : comp)
            for (int k = 0; k < 2; ++k) {
                NodeId n = net.E(e).nd[k];
                if (!net.N(n).hybrid) continue;
                // QUIRK: SNaQ counts a hybrid node once per incident edge of
                // the component (as in the Julia code), not once per component.
                int c = ++counts[n];
                if (c > 2) return false;
            }
    }
    return true;
}

void removehybrid(Network& net, NodeId hyb, bool minor) {
    EdgeId rm = minor ? net.getparentedgeminor(hyb) : net.getparentedge(hyb);
    deletehybridedge(net, rm);
}

void removeexcessreticulations(Network& net, int hmax) {
    while (net.numhybrids > hmax) {
        NodeId best = net.hybrid[0];
        double bg = net.E(net.getparentedgeminor(best)).gamma;
        for (NodeId h : net.hybrid) {
            double g = net.E(net.getparentedgeminor(h)).gamma;
            if (g < bg) {
                bg = g;
                best = h;
            }
        }
        removehybrid(net, best);
        semidirectnetwork(net);
    }
}

}  // namespace csnaq
