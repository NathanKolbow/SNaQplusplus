#include "network.hpp"

#include <algorithm>
#include <unordered_map>

namespace csnaq {

void Edge::remove_node(NodeId n) {
    if (nnode > 0 && nd[0] == n) {
        nd[0] = nd[1];
        nd[1] = NONE;
        --nnode;
    } else if (nnode > 1 && nd[1] == n) {
        nd[1] = NONE;
        --nnode;
    } else {
        throw SnaqError("node " + std::to_string(n) + " not in edge");
    }
}

NodeId Network::new_node(int number, bool isleaf, bool ishybrid) {
    Node n;
    n.number = number;
    n.leaf = isleaf;
    n.hybrid = ishybrid;
    n.booln1 = ishybrid;
    npool.push_back(std::move(n));
    return static_cast<NodeId>(npool.size() - 1);
}

EdgeId Network::new_edge(int number, double length) {
    Edge e;
    e.number = number;
    e.length = length;
    epool.push_back(e);
    return static_cast<EdgeId>(epool.size() - 1);
}

EdgeId Network::new_edge(int number, double length, bool hyb, double gamma, bool ismajor) {
    Edge e;
    e.number = number;
    e.length = length;
    e.hybrid = hyb;
    e.gamma = hyb ? gamma : 1.0;
    e.ismajor = ismajor;
    e.containroot = !hyb;
    epool.push_back(e);
    return static_cast<EdgeId>(epool.size() - 1);
}

NodeList Network::getchildren(NodeId n) const {
    NodeList out;
    for (EdgeId e : N(n).edge)
        if (isparentof(n, e)) out.push_back(getchild(e));
    return out;
}

NodeList Network::getparents(NodeId n) const {
    NodeList out;
    for (EdgeId e : N(n).edge)
        if (ischildof(n, e)) out.push_back(getparent(e));
    return out;
}

NodeId Network::getparent_node(NodeId n) const {
    for (EdgeId e : N(n).edge)
        if (E(e).ismajor && ischildof(n, e)) return getparent(e);
    throw SnaqError("could not find major parent of node " + std::to_string(N(n).number));
}

EdgeId Network::getparentedge(NodeId n) const {
    for (EdgeId e : N(n).edge)
        if (E(e).ismajor && ischildof(n, e)) return e;
    throw SnaqError("node " + std::to_string(N(n).number) + " has no major parent");
}

EdgeId Network::getparentedgeminor(NodeId n) const {
    for (EdgeId e : N(n).edge)
        if (!E(e).ismajor && n == E(e).nd[E(e).ischild1 ? 0 : 1]) return e;
    throw SnaqError("node " + std::to_string(N(n).number) + " has no minor parent");
}

EdgeId Network::getchildedge(NodeId n) const {
    EdgeId found = NONE;
    int count = 0;
    for (EdgeId e : N(n).edge)
        if (isparentof(n, e)) {
            if (count == 0) found = e;
            ++count;
        }
    if (count != 1)
        throw SnaqError("node number " + std::to_string(N(n).number) + " has " + std::to_string(count) +
                        " children instead of 1 child");
    return found;
}

EdgeId Network::getpartneredge(EdgeId e, NodeId n) const {
    for (EdgeId e2 : N(n).edge)
        if (E(e2).hybrid && e2 != e && n == getchild(e2)) return e2;
    throw SnaqError("did not find a partner for edge " + std::to_string(E(e).number));
}

EdgeId Network::getconnectingedge(NodeId a, NodeId b) const {
    for (EdgeId e1 : N(a).edge)
        for (EdgeId e2 : N(b).edge)
            if (e1 == e2) return e1;
    return NONE;
}

bool Network::isexternal(EdgeId e) const {
    const Edge& ed = E(e);
    for (int i = 0; i < ed.nnode; ++i)
        if (N(ed.nd[i]).leaf) return true;
    return false;
}

bool Network::isTree() const {
    if (numhybrids != static_cast<int>(hybrid.size()))
        throw SnaqError("numhybrids does not match to length of net.hybrid");
    return numhybrids == 0;
}

int Network::index_of_node(NodeId n) const {
    for (size_t i = 0; i < node.size(); ++i)
        if (node[i] == n) return static_cast<int>(i);
    return -1;
}

int Network::index_of_edge(EdgeId e) const {
    for (size_t i = 0; i < edge.size(); ++i)
        if (edge[i] == e) return static_cast<int>(i);
    return -1;
}

NodeId Network::find_node_by_number(int number) const {
    for (NodeId n : node)
        if (N(n).number == number) return n;
    return NONE;
}

std::vector<std::string> Network::tiplabels() const {
    std::vector<std::string> out;
    out.reserve(leaf.size());
    for (NodeId l : leaf) out.push_back(N(l).name);
    return out;
}

int Network::max_edge_number() const {
    if (edge.empty()) throw SnaqError("maximum over an empty collection of edges");
    int m = E(edge[0]).number;
    for (EdgeId e : edge) m = std::max(m, E(e).number);
    return m;
}

int Network::max_node_number() const {
    if (node.empty()) throw SnaqError("maximum over an empty collection of nodes");
    int m = N(node[0]).number;
    for (NodeId n : node) m = std::max(m, N(n).number);
    return m;
}

// ------------------------------------------------------------------------

void Network::pushEdge(EdgeId e) {
    edge.push_back(e);
    ++numedges;
}

void Network::pushNode(NodeId n) {
    node.push_back(n);
    ++numnodes;
    if (N(n).leaf) {
        ++numtaxa;
        leaf.push_back(n);
    }
    if (N(n).hybrid) pushHybrid(n);
}

void Network::pushHybrid(NodeId n) {
    if (!N(n).hybrid)
        throw SnaqError("node " + std::to_string(N(n).number) + " is not hybrid, so cannot be pushed in net.hybrid");
    hybrid.push_back(n);
    ++numhybrids;
}

void Network::setNode(EdgeId e, NodeId n) {
    Edge& ed = E(e);
    ed.push_node(n);
    if (ed.nnode == 1) {
        if (ed.hybrid) ed.ischild1 = N(n).hybrid;
        ed.boole1 = !N(n).leaf;
    } else {
        if (N(n).leaf) {
            if (N(ed.nd[0]).leaf) throw SnaqError("edge " + std::to_string(ed.number) + " has two leaves");
            ed.boole1 = false;
        } else if (ed.hybrid) {
            if (N(n).hybrid) {
                ed.ischild1 = false;
            } else {
                if (!N(ed.nd[0]).hybrid)
                    throw SnaqError("hybrid edge " + std::to_string(ed.number) + " has no hybrid nodes");
                ed.ischild1 = true;
            }
        }
    }
}

void Network::setNode(EdgeId e, NodeId n1, NodeId n2) {
    Edge& ed = E(e);
    ed.nd[0] = n1;
    ed.nd[1] = n2;
    ed.nnode = 2;
    if (ed.hybrid) {
        if (N(n1).hybrid) {
            ed.ischild1 = true;
        } else {
            if (!N(n2).hybrid) throw SnaqError("hybrid edge without hybrid node");
            ed.ischild1 = false;
        }
    }
    ed.boole1 = !(N(n1).leaf || N(n2).leaf);
}

void Network::setEdge(NodeId n, EdgeId e) {
    Node& nd = N(n);
    nd.edge.push_back(e);
    bool any = false;
    for (EdgeId x : nd.edge) any |= E(x).hybrid;
    nd.booln1 = any;
}

void Network::removeEdge(NodeId n, EdgeId e) {
    Node& nd = N(n);
    int i = nd.edge.find(e);
    if (i < 0)
        throw SnaqError("edge " + std::to_string(E(e).number) + " not in node " + std::to_string(nd.number));
    nd.edge.erase_at(static_cast<size_t>(i));
    bool any = false;
    for (EdgeId x : nd.edge) any |= E(x).hybrid;
    nd.booln1 = any;
}

void Network::removeNode(NodeId n, EdgeId e) { E(e).remove_node(n); }

void Network::deleteEdge(EdgeId e) {
    int i = index_of_edge(e);
    if (i < 0) throw SnaqError("edge " + std::to_string(E(e).number) + " not in network: can't delete");
    edge.erase(edge.begin() + i);
    --numedges;
}

void Network::deleteNode(NodeId n) {
    int index = index_of_node(n);
    if (index < 0) throw SnaqError("Node " + std::to_string(N(n).number) + " not in network");
    node.erase(node.begin() + index);
    --numnodes;
    if (rooti == index)
        rooti = 0;  // arbitrary, as in PhyloNetworks
    else if (rooti > index)
        --rooti;
    if (N(n).hybrid) removeHybrid(n);
    if (N(n).leaf) removeLeaf(n);
}

void Network::removeHybrid(NodeId n) {
    if (!N(n).hybrid)
        throw SnaqError("cannot delete node " + std::to_string(N(n).number) +
                        " from net.hybrid because it is not hybrid");
    auto it = std::find(hybrid.begin(), hybrid.end(), n);
    if (it == hybrid.end())
        throw SnaqError("hybrid node " + std::to_string(N(n).number) + " not in the network's list of hybrids");
    hybrid.erase(it);
    --numhybrids;
}

void Network::removeLeaf(NodeId n) {
    if (!N(n).leaf)
        throw SnaqError("cannot delete node " + std::to_string(N(n).number) + " from net.leaf because it is not leaf");
    auto it = std::find(leaf.begin(), leaf.end(), n);
    if (it == leaf.end()) throw SnaqError("leaf node " + std::to_string(N(n).number) + " not in network");
    leaf.erase(it);
    --numtaxa;
}

Network Network::deepcopynetwork() const {
    Network c;
    c.npool.reserve(node.size());
    c.epool.reserve(edge.size());
    std::vector<NodeId> nmap(npool.size(), NONE);
    // QUIRK: SNaQ.jl keys its node map by node *number*; numbers are unique
    // whenever a copy is taken, so mapping by identity is equivalent.
    for (NodeId n : node) {
        const Node& src = N(n);
        NodeId id = c.new_node(src.number, src.leaf, src.hybrid);
        c.N(id).name = src.name;
        nmap[static_cast<size_t>(n)] = id;
        c.node.push_back(id);
    }
    for (EdgeId e : edge) {
        const Edge& src = E(e);
        NodeId n1 = nmap[static_cast<size_t>(src.nd[0])];
        NodeId n2 = nmap[static_cast<size_t>(src.nd[1])];
        // Edge(number, length, hybrid, gamma, [n1, n2])
        Edge ne;
        ne.number = src.number;
        ne.length = src.length;
        ne.hybrid = src.hybrid;
        ne.gamma = src.hybrid ? src.gamma : 1.0;  // QUIRK: tree edges get γ = 1 in copies
        ne.nd[0] = n1;
        ne.nd[1] = n2;
        ne.nnode = 2;
        ne.ischild1 = src.ischild1;
        ne.ismajor = src.ismajor;
        ne.containroot = src.containroot;
        c.epool.push_back(ne);
        EdgeId id = static_cast<EdgeId>(c.epool.size() - 1);
        c.N(n1).edge.push_back(id);
        c.N(n2).edge.push_back(id);
        c.edge.push_back(id);
    }
    for (NodeId h : hybrid) c.hybrid.push_back(nmap[static_cast<size_t>(h)]);
    for (NodeId l : leaf) c.leaf.push_back(nmap[static_cast<size_t>(l)]);
    c.numtaxa = numtaxa;
    c.numnodes = numnodes;
    c.numedges = numedges;
    c.numhybrids = numhybrids;
    c.rooti = rooti;
    c.names = names;
    c.isrooted = isrooted;
    c.fscore = fscore;
    return c;
}

}  // namespace csnaq
