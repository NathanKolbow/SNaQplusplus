#include "compare.hpp"

#include <algorithm>
#include <limits>
#include <set>
#include <string>
#include <vector>

#include "network_manip.hpp"

namespace snaqpp {

// ------------------------------------------------------------------------
// hardwired clusters (PhyloNetworks compareNetworks.jl)
// ------------------------------------------------------------------------

namespace {

using Row = std::vector<int>;  // [edge number, cluster indicators..., (type)]

int find_taxon(const std::vector<std::string>& S, const std::string& name) {
    auto it = std::find(S.begin(), S.end(), name);
    if (it == S.end()) throw SnaqError("leaf " + name + " not in taxon list");
    return static_cast<int>(it - S.begin());
}

void traverseTree2Matrix(const Network& T, NodeId node, EdgeId edge, std::vector<Row>& M, size_t& ie,
                         const std::vector<std::string>& S) {
    NodeId child = T.getOtherNode(edge, node);
    size_t indedge = ie;
    M[indedge][0] = T.E(edge).number;
    ++ie;
    for (EdgeId e : T.N(child).edge) {
        if (e == edge) continue;
        NodeId gc = T.getOtherNode(e, child);
        if (T.N(gc).leaf) {
            M[indedge][static_cast<size_t>(find_taxon(S, T.N(gc).name)) + 1] = 1;
        } else {
            size_t inde = ie;
            traverseTree2Matrix(T, child, e, M, ie, S);
            for (size_t k = 1; k < M[indedge].size(); ++k) M[indedge][k] |= M[inde][k];
        }
    }
}

std::vector<Row> tree2Matrix(const Network& T, const std::vector<std::string>& S, bool rooted) {
    size_t ne = T.edge.size() - static_cast<size_t>(T.numtaxa);
    NodeId root = T.getroot();
    if (T.N(root).leaf) ++ne;
    std::vector<Row> M(ne, Row(S.size() + 1, 0));
    size_t ie = 0;
    for (EdgeId e : T.N(root).edge)
        if (!T.N(T.getOtherNode(e, root)).leaf) traverseTree2Matrix(T, root, e, M, ie, S);
    if (!rooted && T.N(root).edge.size() < 3 && !M.empty()) M.erase(M.begin());
    return M;
}

int hardwiredclusters_rec(const Network& net, NodeId node, EdgeId edge, std::vector<Row>& M, size_t& ie,
                          const std::vector<std::string>& S) {
    NodeId child = net.getOtherNode(edge, node);
    if (net.N(child).leaf) {
        find_taxon(S, net.N(child).name);
        return -1;
    }
    EdgeId partner = NONE;
    if (net.E(edge).hybrid) {
        if (net.getchild(edge) != child) throw SnaqError("inconsistency during network traversal");
        bool partnerVisited = true;
        int indpartner = -1;
        for (EdgeId e : net.N(child).edge)
            if (net.E(e).hybrid && e != edge && net.getchild(e) == child) {
                partner = e;
                for (size_t r = 0; r < ie; ++r)
                    if (M[r][0] == net.E(e).number) {
                        indpartner = static_cast<int>(r);
                        break;
                    }
                if (indpartner < 0) partnerVisited = false;
                break;
            }
        if (partner == NONE) throw SnaqError("partner hybrid edge not found");
        if (partnerVisited) return indpartner;
    }
    size_t indedge = ie;
    M[indedge][0] = net.E(edge).number;
    M[indedge].back() = net.E(edge).hybrid ? 11 : 10;
    ++ie;
    for (EdgeId e : net.N(child).edge) {
        if (e == edge || (net.E(edge).hybrid && e == partner)) continue;
        NodeId gc = net.getOtherNode(e, child);
        if (net.N(gc).leaf) {
            M[indedge][static_cast<size_t>(find_taxon(S, net.N(gc).name)) + 1] = 1;
        } else {
            int inde = hardwiredclusters_rec(net, child, e, M, ie, S);
            for (size_t k = 1; k + 1 < M[indedge].size(); ++k) M[indedge][k] |= M[static_cast<size_t>(inde)][k];
        }
    }
    return static_cast<int>(indedge);
}

std::vector<Row> hardwiredclusters(const Network& net, const std::vector<std::string>& S) {
    size_t ne = net.edge.size() - static_cast<size_t>(net.numtaxa) - net.hybrid.size();
    NodeId root = net.getroot();
    if (net.N(root).leaf) ++ne;
    std::vector<Row> M(ne, Row(S.size() + 2, 0));
    size_t ie = 0;
    for (EdgeId e : net.N(root).edge) hardwiredclusters_rec(net, root, e, M, ie, S);
    return M;
}

int hwcd_rooted_or_trees(const Network& net1, const Network& net2, bool rooted) {
    bool bothtrees = net1.numhybrids == 0 && net2.numhybrids == 0;
    std::vector<std::string> taxa = net1.tiplabels();
    std::sort(taxa.begin(), taxa.end());
    auto t2 = net2.tiplabels();
    std::set<std::string> s2(t2.begin(), t2.end());
    for (const auto& t : taxa)
        if (!s2.count(t)) throw SnaqError("net1 and net2 do not share the same taxon set.");
    std::vector<Row> M1, M2;
    if (bothtrees) {
        M1 = tree2Matrix(net1, taxa, rooted);
        M2 = tree2Matrix(net2, taxa, rooted);
    } else {
        M1 = hardwiredclusters(net1, taxa);
        M2 = hardwiredclusters(net2, taxa);
    }
    int dis = 0;
    std::vector<size_t> n2ci(M2.size());
    for (size_t k = 0; k < M2.size(); ++k) n2ci[k] = k;
    for (const Row& r1 : M1) {
        Row a(r1.begin() + 1, r1.end());
        Row m1(a.size());
        for (size_t k = 0; k < a.size(); ++k) m1[k] = 1 - a[k];
        bool found = false;
        for (size_t j = n2ci.size(); j-- > 0;) {
            const Row& r2 = M2[n2ci[j]];
            Row b(r2.begin() + 1, r2.end());
            if (a == b || (!rooted && m1 == b)) {
                found = true;
                n2ci.erase(n2ci.begin() + static_cast<long>(j));
                break;
            }
        }
        if (!found) ++dis;
    }
    return dis + dis + static_cast<int>(M2.size()) - static_cast<int>(M1.size());
}

}  // namespace

int hardwiredclusterdistance(const Network& n1, const Network& n2, bool rooted) {
    bool bothtrees = n1.numhybrids == 0 && n2.numhybrids == 0;
    if (rooted || bothtrees) return hwcd_rooted_or_trees(n1, n2, rooted);
    Network net1 = n1.clone(), net2 = n2.clone();
    removedegree2nodes(net1);
    removedegree2nodes(net2);
    auto roots = [](Network& net) {
        std::vector<int> r;
        for (NodeId n : net.node)
            if (!net.N(n).leaf) r.push_back(net.N(n).number);
        for (size_t i = r.size(); i-- > 0;) {
            try {
                rootatnode_number(net, r[i]);
            } catch (const RootMismatch&) {
                r.erase(r.begin() + static_cast<long>(i));
            }
        }
        return r;
    };
    std::vector<int> r1 = roots(net1), r2 = roots(net2);
    int best = std::numeric_limits<int>::max();
    for (int a : r1) {
        rootatnode_number(net1, a);
        for (int b : r2) {
            rootatnode_number(net2, b);
            best = std::min(best, hwcd_rooted_or_trees(net1, net2, true));
        }
    }
    return best;
}

}  // namespace snaqpp
