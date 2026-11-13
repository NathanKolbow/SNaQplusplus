#include "network_manip.hpp"

#include <algorithm>
#include <regex>
#include <set>

namespace csnaq {

// ------------------------------------------------------------------------
// directedges! and helpers (manipulateNet.jl)
// ------------------------------------------------------------------------

static void traverseDirectEdges(Network& net, NodeId node, EdgeId edge, bool containroot) {
    Edge& ed = net.E(edge);
    if (ed.hybrid && node == net.getchild(edge))
        throw RootMismatch("direction (ischild1) of hybrid edge " + std::to_string(ed.number) +
                           " conflicts with the root.\nischild1 and containroot were updated for a "
                           "subset of edges in the network only.");
    NodeId cn;
    if (node == ed.nd[0]) {
        ed.ischild1 = false;
        cn = ed.nd[1];
    } else {
        ed.ischild1 = true;
        cn = ed.nd[0];
    }
    ed.containroot = containroot;
    if (!net.N(cn).leaf && (!ed.hybrid || ed.ismajor)) {
        if (ed.hybrid) containroot = false;
        int nchildren = 0;
        // iterate over a copy: recursion does not modify cn's edge list,
        // but ischild1 of its edges is read at iteration time, as in Julia.
        for (size_t k = 0; k < net.N(cn).edge.size(); ++k) {
            EdgeId e = net.N(cn).edge[k];
            if (e == edge) continue;
            if (net.E(e).hybrid && cn == net.getchild(e)) continue;
            traverseDirectEdges(net, cn, e, containroot);
            ++nchildren;
        }
        if (nchildren == 0)
            throw RootMismatch("non-leaf node " + std::to_string(net.N(cn).number) +
                               " had 0 children.\nCould be a hybrid whose parents' direction conflicts with "
                               "the root.\nischild1 and containroot were updated for a subset of edges in the "
                               "network only.");
    }
}

void directedges(Network& net, bool checkMajor) {
    if (checkMajor) {
        for (NodeId n : net.node) {
            int nparents = 0, nmajor = 0;
            for (EdgeId e : net.N(n).edge) {
                if (net.E(e).hybrid && n == net.getchild(e)) {
                    ++nparents;
                    if (net.E(e).ismajor) ++nmajor;
                }
            }
            if (nparents == 1)
                throw SnaqError("node " + std::to_string(net.N(n).number) + " has exactly 1 hybrid parent edge");
            if (!(nparents == 0 || nmajor == 1))
                throw SnaqError("hybrid node " + std::to_string(net.N(n).number) + " has 0 or 2+ major hybrid parents");
        }
    }
    NodeId root = net.getroot();
    for (size_t k = 0; k < net.N(root).edge.size(); ++k) traverseDirectEdges(net, root, net.N(root).edge[k], true);
    net.isrooted = true;
}

void preorder(Network& net) {
    if (!net.isrooted) throw SnaqError("net needs to be rooted for preorder!, run root functions or directedges!");
    net.vec_node.clear();
    std::vector<char> visited(net.node.size(), 0);
    std::vector<int> queue{net.rooti};
    while (!queue.empty()) {
        int currind = queue.back();
        queue.pop_back();
        NodeId curr = net.node[static_cast<size_t>(currind)];
        if (visited[static_cast<size_t>(currind)]) continue;
        visited[static_cast<size_t>(currind)] = 1;
        net.vec_node.push_back(curr);
        for (EdgeId e : net.N(curr).edge) {
            if (net.getparent(e) != curr) continue;
            NodeId cn = net.getchild(e);
            if (!net.E(e).hybrid) {
                queue.push_back(net.index_of_node(cn));
            } else {
                EdgeId e2 = net.getpartneredge(e, cn);
                NodeId parent = net.getparent(e2);
                if (visited[static_cast<size_t>(net.index_of_node(parent))]) queue.push_back(net.index_of_node(cn));
            }
        }
    }
    for (char v : visited)
        if (!v)
            throw RootMismatch("Not all nodes were visited after 'preorder!'.\nWere edges directed away from the "
                               "root? are there multiple roots?");
}

void allowrootbelow(Network& net, EdgeId e) {
    net.E(e).containroot = true;
    if (net.E(e).hybrid) return;
    allowrootbelow(net, net.getchild(e), e);
}

void allowrootbelow(Network& net, NodeId n, EdgeId pe) {
    for (size_t k = 0; k < net.N(n).edge.size(); ++k) {
        EdgeId ce = net.N(n).edge[k];
        if (ce == pe) continue;
        if (net.getparent(ce) != n) continue;
        allowrootbelow(net, ce);
    }
}

void allowrootbelow(Network& net) {
    NodeId rn = net.getroot();
    for (size_t k = 0; k < net.N(rn).edge.size(); ++k) {
        EdgeId e = net.N(rn).edge[k];
        if (net.E(e).containroot) allowrootbelow(net, net.getchild(e), e);
    }
}

void norootbelow(Network& net, EdgeId e) {
    if (!net.E(e).containroot) return;
    net.E(e).containroot = false;
    NodeId cn = net.getchild(e);
    for (size_t k = 0; k < net.N(cn).edge.size(); ++k) {
        EdgeId ce = net.N(cn).edge[k];
        if (ce == e) continue;
        if (net.getparent(ce) != cn) continue;
        norootbelow(net, ce);
    }
}

static bool isdescendant_rec(const Network& net, std::vector<int>& visited, NodeId des, EdgeId e) {
    NodeId n = net.getchild(e);
    if (n == des) return true;
    if (net.N(n).hybrid) {
        int num = net.N(n).number;
        if (std::find(visited.begin(), visited.end(), num) != visited.end()) return false;
        visited.push_back(num);
    }
    for (EdgeId ce : net.N(n).edge)
        if (net.isparentof(n, ce) && isdescendant_rec(net, visited, des, ce)) return true;
    return false;
}

bool isdescendant(const Network& net, NodeId des, NodeId anc) {
    std::vector<int> visited;
    for (EdgeId e : net.N(anc).edge) {
        if (anc == net.getchild(e)) continue;
        if (isdescendant_rec(net, visited, des, e)) return true;
    }
    return false;
}

static bool isdescendant_undirected_rec(const Network& net, std::vector<int>& visited, NodeId des, NodeId anc,
                                        EdgeId parentedge) {
    for (EdgeId e : net.N(anc).edge) {
        if (e == parentedge) continue;
        if (net.E(e).hybrid && net.getparent(e) != anc) continue;
        NodeId n = net.getOtherNode(e, anc);
        if (n == des) return true;
        if (net.N(n).hybrid) {
            int num = net.N(n).number;
            if (std::find(visited.begin(), visited.end(), num) != visited.end()) continue;
            visited.push_back(num);
        }
        if (isdescendant_undirected_rec(net, visited, des, n, e)) return true;
    }
    return false;
}

bool isdescendant_undirected(const Network& net, NodeId des, NodeId anc, EdgeId parentedge) {
    std::vector<int> visited;
    return isdescendant_undirected_rec(net, visited, des, anc, parentedge);
}

// ------------------------------------------------------------------------
// deletion
// ------------------------------------------------------------------------

void deleteleaf_node(Network& net, NodeId n, const DeleteOpts& o) {
    if (!net.N(n).leaf) throw SnaqError("node number " + std::to_string(net.N(n).number) + " is not a leaf.");
    deleteleaf_number(net, net.N(n).number, o);
}

void deleteleaf_name(Network& net, const std::string& name, const DeleteOpts& o) {
    int found = -1, count = 0;
    for (size_t i = 0; i < net.node.size(); ++i)
        if (net.N(net.node[i]).name == name) {
            if (count == 0) found = static_cast<int>(i);
            ++count;
        }
    if (count == 0) throw SnaqError("node named " + name + " was not found in the network.");
    if (count > 1) throw SnaqError("several nodes were found with name " + name + ".");
    deleteleaf_index(net, found, o);
}

void deleteleaf_number(Network& net, int nodeNumber, const DeleteOpts& o) {
    int i = -1;
    for (size_t k = 0; k < net.node.size(); ++k)
        if (net.N(net.node[k]).number == nodeNumber) {
            i = static_cast<int>(k);
            break;
        }
    if (i < 0)
        throw SnaqError("cannot delete leaf number " + std::to_string(nodeNumber) + " because it is not part of net");
    deleteleaf_index(net, i, o);
}

void deleteleaf_index(Network& net, int i, const DeleteOpts& o) {
    if (i >= static_cast<int>(net.node.size()))
        throw SnaqError("node index too large: the network only has " + std::to_string(net.node.size()) + " nodes.");
    NodeId nodei = net.node[static_cast<size_t>(i)];
    size_t nodeidegree = net.N(nodei).edge.size();
    if (nodeidegree == 0) {
        if (net.node.size() != 1) throw SnaqError("leaf has no edge but network has more than 1 node");
        warn("Only 1 node. Removing it: the network will be empty");
        net.deleteNode(nodei);
        return;
    } else if (nodeidegree == 1) {
        EdgeId pe = net.N(nodei).edge[0];
        NodeId pn = net.getOtherNode(pe, nodei);
        if (net.rooti == i && o.keeporiginalroot) return;
        if (net.N(pn).leaf) {
            net.rooti = i;
            if (net.edge.size() != 1)
                throw SnaqError("neighbor of degree-1 node " + net.N(nodei).name + " is a leaf, but network has " +
                                std::to_string(net.edge.size()) + " edges (instead of 1).");
            if (net.N(pn).edge.size() != 1) throw SnaqError("neighbor of degree-1 node is a leaf with 2+ edges");
            return;
        }
        net.removeNode(pn, pe);
        net.removeEdge(pn, pe);
        net.deleteEdge(pe);
        if (net.rooti == i) net.rooti = net.index_of_node(pn);
        net.deleteNode(nodei);
        deleteleaf_number(net, net.N(pn).number, o);
        return;
    } else if (nodeidegree > 2) {
        return;
    }
    // degree 2
    if (i == net.rooti && (o.keeporiginalroot || !o.unroot)) return;
    EdgeId e1 = net.N(nodei).edge[0];
    EdgeId e2 = net.N(nodei).edge[1];
    if (net.E(e1).hybrid && net.E(e2).hybrid) {
        NodeId cn = net.getchild(e1);
        NodeId cn2 = net.getchild(e2);
        if (!(nodei == cn && nodei == cn2)) {
            if (!(i == net.rooti && cn == cn2))
                throw SnaqError("after removing descendants, node " + std::to_string(net.N(nodei).number) +
                                " has 2 hybrid edges but is not the child of both.");
            if (!net.N(cn).hybrid) throw SnaqError("child node of hybrid edges should be a hybrid.");
            for (EdgeId e : net.N(cn).edge)
                if (net.getchild(e) == cn && e != e1 && e != e2)
                    throw SnaqError("root has 2 hybrid edges, but their common child has an extra parent");
            net.removeEdge(cn, e1);
            net.removeEdge(cn, e2);
            net.removeHybrid(cn);
            net.N(cn).hybrid = false;
            net.E(e1).nnode = 0;
            net.E(e2).nnode = 0;
            net.deleteEdge(e1);
            net.deleteEdge(e2);
            net.N(nodei).edge.clear();
            net.deleteNode(nodei);
            net.rooti = net.index_of_node(cn);
            deleteleaf_index(net, net.rooti, o);
            return;
        }
        NodeId p1 = net.getparent(e1);
        NodeId p2 = net.getparent(e2);
        bool sameparent = (p1 == p2);
        net.removeNode(p1, e1);
        net.removeNode(p2, e2);
        net.removeEdge(p1, e1);
        net.removeEdge(p2, e2);
        net.deleteEdge(e1);
        net.deleteEdge(e2);
        if (net.rooti == i) net.rooti = net.index_of_node(p1);
        net.deleteNode(nodei);
        deleteleaf_number(net, net.N(p1).number, o);
        if (!sameparent) {
            int p2num = net.N(p2).number;
            int p2idx = -1;
            for (size_t k = 0; k < net.node.size(); ++k)
                if (net.N(net.node[k]).number == p2num) {
                    p2idx = static_cast<int>(k);
                    break;
                }
            if (p2idx >= 0) deleteleaf_index(net, p2idx, o);
        }
    } else if (!o.nofuse) {
        EdgeId fe = fuseedgesat(net, i, o.multgammas);
        if (o.simplify && net.E(fe).hybrid) {
            NodeId cn = net.getchild(fe);
            EdgeId fe2 = net.getpartneredge(fe, cn);
            NodeId pn = net.getparent(fe);
            if (pn == net.getparent(fe2)) {
                net.E(fe2).hybrid = false;
                net.E(fe2).ismajor = true;
                net.E(fe2).gamma = addBL(net.E(fe).gamma, net.E(fe2).gamma);
                net.removeEdge(pn, fe);
                net.removeEdge(cn, fe);
                net.deleteEdge(fe);
                net.removeHybrid(cn);
                net.N(cn).hybrid = false;
                deleteleaf_number(net, net.N(cn).number, o);
                int pnnum = net.N(pn).number;
                int pnidx = -1;
                for (size_t k = 0; k < net.node.size(); ++k)
                    if (net.N(net.node[k]).number == pnnum) {
                        pnidx = static_cast<int>(k);
                        break;
                    }
                if (pnidx >= 0) deleteleaf_index(net, pnidx, o);
            }
        }
    }
}

void deletehybridedge(Network& net, EdgeId edge, bool nofuse, bool unroot, bool multgammas, bool simplify,
                      bool keeporiginalroot) {
    if (!net.E(edge).hybrid)
        throw SnaqError("edge " + std::to_string(net.E(edge).number) + " has to be hybrid for deletehybridedge!");
    NodeId n1 = net.getchild(edge);
    if (!net.N(n1).hybrid) throw SnaqError("child node of hybrid edge should be a hybrid.");
    size_t n1degree = net.N(n1).edge.size();
    NodeId n2 = net.getparent(edge);
    size_t n2degree = net.N(n2).edge.size();
    bool delete_n1_recursively = false;
    if (n1degree < 3) {
        throw SnaqError("node " + std::to_string(net.N(n1).number) + " has " + std::to_string(n1degree) +
                        " edges instead of 3+");
    } else if (n1degree == 3 && !nofuse) {
        delete_n1_recursively = true;
        EdgeId pe = NONE, ce = NONE;
        for (EdgeId e : net.N(n1).edge) {
            if (net.E(e).hybrid && e != edge && n1 == net.getchild(e)) pe = e;
            if (!net.E(e).hybrid || n1 == net.getparent(e)) ce = e;
        }
        NodeId pn = net.getparent(pe);
        bool atRoot = (net.getroot() == n1);
        if (net.E(pe).containroot) allowrootbelow(net, ce);
        net.E(ce).length = addBL(net.E(ce).length, net.E(pe).length);
        if (multgammas) net.E(ce).gamma = multiplygammas(net.E(ce).gamma, net.E(pe).gamma);
        net.removeNode(n1, ce);
        net.setNode(ce, pn);
        net.E(ce).ischild1 = true;
        net.setEdge(pn, ce);
        net.removeEdge(pn, pe);
        net.removeEdge(n1, pe);
        net.removeEdge(n1, ce);
        net.deleteEdge(pe);
        net.removeHybrid(n1);
        net.N(n1).hybrid = false;
        net.E(edge).hybrid = false;
        net.E(edge).ismajor = true;
        if (atRoot) {
            int i = net.index_of_node(pn);
            if (i < 0) throw SnaqError("node not in net!");
            net.rooti = i;
        }
    } else {
        EdgeId pe = net.getpartneredge(edge, n1);
        if (!net.E(pe).ismajor) net.E(pe).ismajor = true;
        net.E(pe).hybrid = false;
        net.removeEdge(n1, edge);
        net.removeHybrid(n1);
        net.N(n1).hybrid = false;
        if (net.E(pe).containroot) allowrootbelow(net, pe);
    }

    int formernumhyb = net.numhybrids;
    DeleteOpts o{nofuse, simplify, unroot, multgammas, keeporiginalroot};
    if (delete_n1_recursively) {
        deleteleaf_number(net, net.N(n1).number, o);
    } else if (n2degree == 1) {
        throw SnaqError("node (parent of hybrid edge to be deleted) has 1 edge only!");
    } else {
        net.removeEdge(n2, edge);
        net.deleteEdge(edge);
        deleteleaf_number(net, net.N(n2).number, o);
    }
    if (net.numhybrids != formernumhyb) allowrootbelow(net);
}

// ------------------------------------------------------------------------
// addition / restructuring
// ------------------------------------------------------------------------

NodeId addleaf(Network& net, NodeId speciesnode, const std::string& leafname, double edgelength) {
    EdgeId ext = net.new_edge(net.max_edge_number() + 1, edgelength);
    net.pushEdge(ext);
    net.setEdge(speciesnode, ext);
    if (net.N(speciesnode).hybrid ||
        (!net.isrootof(speciesnode) && !net.E(net.getparentedge(speciesnode)).containroot))
        net.E(ext).containroot = false;
    NodeId newleaf = net.new_node(net.max_node_number() + 1, true, false);
    net.N(newleaf).edge.push_back(ext);
    net.N(newleaf).booln1 = net.E(ext).hybrid;
    net.N(newleaf).name = leafname;
    net.setNode(ext, newleaf, speciesnode);
    if (net.N(speciesnode).leaf) {
        auto it = std::find(net.leaf.begin(), net.leaf.end(), speciesnode);
        net.leaf.erase(it);
        net.N(speciesnode).leaf = false;
        --net.numtaxa;
    }
    net.pushNode(newleaf);
    return newleaf;
}

NodeId addleaf_on_edge(Network& net, EdgeId startingedge, const std::string& leafname, double edgelength) {
    NodeId newnode = breakedge(net, startingedge).first;
    return addleaf(net, newnode, leafname, edgelength);
}

std::pair<NodeId, EdgeId> breakedge(Network& net, EdgeId edge) {
    NodeId pn = net.getparent(edge);
    net.removeEdge(pn, edge);
    net.removeNode(pn, edge);
    int max_edge = net.max_edge_number();
    int max_node = net.max_node_number();
    EdgeId newedge = net.new_edge(max_edge + 1, -1.0);
    {
        Edge& ne = net.E(newedge);
        ne.hybrid = false;
        ne.y = -1.0;
        ne.z = -1.0;
        ne.gamma = 1.0;
        ne.ischild1 = true;
        ne.ismajor = true;
        ne.inte1 = -1;
        ne.containroot = net.E(edge).containroot;
        ne.boole1 = true;
        ne.boole2 = false;
    }
    NodeId newnode = net.new_node(max_node + 1, false, false);
    net.N(newnode).edge.push_back(edge);
    net.N(newnode).edge.push_back(newedge);
    net.N(newnode).booln1 = net.E(edge).hybrid;
    net.setNode(edge, newnode);
    net.E(edge).ischild1 = true;
    net.setNode(newedge, newnode);
    net.setEdge(pn, newedge);
    net.setNode(newedge, pn);
    double oldlen = net.E(edge).length;
    if (oldlen != -1.0) {
        double len_parent = 0.5 * oldlen;
        net.E(newedge).length = len_parent;
        net.E(edge).length = oldlen - len_parent;
    }
    net.pushEdge(newedge);
    net.pushNode(newnode);
    return {newnode, newedge};
}

EdgeId fuseedgesat(Network& net, int i, bool multgammas) {
    if (i >= static_cast<int>(net.node.size())) throw SnaqError("node index too large in fuseedgesat!");
    NodeId nodei = net.node[static_cast<size_t>(i)];
    if (net.N(nodei).edge.size() != 2)
        throw SnaqError("can't fuse edges at node number " + std::to_string(net.N(nodei).number) +
                        ": connected to " + std::to_string(net.N(nodei).edge.size()) + " edges.");
    EdgeId ea = net.N(nodei).edge[0], eb = net.N(nodei).edge[1];
    if (net.E(ea).hybrid && net.E(eb).hybrid)
        throw SnaqError("can't fuse edges at node number " + std::to_string(net.N(nodei).number) +
                        ": connected to exactly 2 hybrid edges");
    int j = (net.E(eb).number > net.E(ea).number) ? 1 : 0;  // argmax: first maximum
    EdgeId pe = net.N(nodei).edge[static_cast<size_t>(j)];
    EdgeId ce = net.N(nodei).edge[static_cast<size_t>(j == 0 ? 1 : 0)];
    if (net.E(pe).hybrid) std::swap(ce, pe);
    bool isnodeiparent = (nodei == net.getparent(ce));
    if (!(!net.E(ce).hybrid || isnodeiparent))
        throw SnaqError("node " + std::to_string(net.N(nodei).number) +
                        " has 1 tree edge and 1 hybrid edge, but is child of the hybrid edge.");
    NodeId pn = net.getOtherNode(pe, nodei);
    net.removeEdge(nodei, ce);
    net.removeNode(nodei, ce);
    net.removeEdge(pn, pe);
    net.removeNode(pn, pe);
    net.setEdge(pn, ce);
    net.setNode(ce, pn);
    net.E(ce).ischild1 = isnodeiparent;
    net.E(ce).length = addBL(net.E(ce).length, net.E(pe).length);
    if (multgammas) net.E(ce).gamma = multiplygammas(net.E(ce).gamma, net.E(pe).gamma);
    if (net.rooti == i) {
        NodeId newroot = pn;
        if (net.N(newroot).leaf && !net.E(ce).hybrid) {
            newroot = net.E(ce).nd[0];
            net.E(ce).ischild1 = false;
        }
        net.rooti = net.index_of_node(newroot);
    }
    net.deleteNode(nodei);
    net.deleteEdge(pe);
    return ce;
}

void removedegree2nodes(Network& net, bool keeproot) {
    NodeId rootnode = net.getroot();
    auto rootin2cycle = [&](NodeId nn) {
        if (!net.isrootof(nn)) return false;
        for (EdgeId e : net.N(nn).edge)
            if (!net.E(e).hybrid) return false;
        return true;
    };
    auto toberemoved = [&](NodeId nn) {
        return keeproot ? (net.N(nn).edge.size() == 2 && nn != rootnode)
                        : (net.N(nn).edge.size() == 2 && !rootin2cycle(nn));
    };
    int ndegree2nodes = 0;
    for (NodeId n : net.node) ndegree2nodes += toberemoved(n) ? 1 : 0;
    for (int it = 0; it < ndegree2nodes; ++it) {
        int found = -1;
        for (size_t k = 0; k < net.node.size(); ++k)
            if (toberemoved(net.node[k])) {
                found = static_cast<int>(k);
                break;
            }
        if (found >= 0) fuseedgesat(net, found);
    }
}

static void shrink2cycleat(Network& net, EdgeId minor, EdgeId major, bool unroot) {
    double g = net.E(minor).gamma;
    if (g == -1.0) g = 0.5;
    net.E(major).length =
        addBL(multiplygammas(g, net.E(minor).length), multiplygammas(1.0 - g, net.E(major).length));
    deletehybridedge(net, minor, false, unroot, false, false, false);
}

static bool shrink3cycleat(Network& net, NodeId /*hybrid*/, EdgeId edge1, EdgeId edge2, NodeId node1, NodeId node2,
                           bool unroot) {
    EdgeId edge3 = NONE;
    for (EdgeId e : net.N(node1).edge) {
        if (e == edge1) continue;
        if (net.getOtherNode(e, node1) == node2) {
            edge3 = e;
            break;
        }
    }
    if (edge3 == NONE) return false;
    if (net.E(edge3).hybrid) {
        if (net.getchild(edge3) == node1) {
            std::swap(node1, node2);
            std::swap(edge1, edge2);
        }
        double g1 = net.E(edge1).gamma;
        double g2g3 = multiplygammas(net.E(edge2).gamma, net.E(edge3).gamma);
        double g1tilde = addBL(g1, g2g3);
        if (g1tilde != -1.0) {
            net.E(edge1).gamma = g1tilde;
            net.E(edge2).gamma = 1.0 - g1tilde;
            net.E(edge1).ismajor = g1tilde >= 0.5;
            net.E(edge2).ismajor = !net.E(edge1).ismajor;
            if (net.E(edge1).length != -1.0 && net.E(edge2).length != -1.0 && net.E(edge3).length != -1.0)
                net.E(edge1).length =
                    (net.E(edge1).length * g1 + (net.E(edge3).length + net.E(edge2).length) * g2g3) / g1tilde;
        }
        deletehybridedge(net, edge3, false, unroot, false, false);
    } else {
        EdgeId edgeB = NONE;
        for (EdgeId e : net.N(node2).edge)
            if (e != edge1 && e != edge3) {
                edgeB = e;
                break;
            }
        double g1t1 = multiplygammas(net.E(edge1).gamma, net.E(edge1).length);
        double t3 = net.E(edge3).length;
        if (g1t1 != -1.0 && net.E(edge2).length != -1.0 && t3 != -1.0) {
            net.E(edgeB).length = addBL(net.E(edgeB).length, net.E(edge1).gamma * t3);
            net.E(edge3).length = t3 * net.E(edge2).gamma;
            net.E(edge2).length = g1t1 + net.E(edge2).gamma * net.E(edge2).length;
        }
        deletehybridedge(net, edge1, false, unroot, false, false);
    }
    return true;
}

bool shrink3cycles(Network& net, bool unroot) {
    bool foundcycle = false;
    int nh = static_cast<int>(net.hybrid.size());
    int ih = nh;  // 1-based as in Julia
    while (ih > 0) {
        NodeId h = net.hybrid[static_cast<size_t>(ih - 1)];
        EdgeId minor = net.getparentedgeminor(h);
        EdgeId major = net.getparentedge(h);
        NodeId pmin = net.getparent(minor);
        NodeId pmaj = net.getparent(major);
        if (pmin == pmaj) {
            foundcycle = true;
            shrink2cycleat(net, minor, major, unroot);
            nh = static_cast<int>(net.hybrid.size());
            ih = nh + 1;
        } else {
            if (shrink3cycleat(net, h, minor, major, pmin, pmaj, unroot)) {
                foundcycle = true;
                nh = static_cast<int>(net.hybrid.size());
                ih = nh + 1;
            }
        }
        ih -= 1;
    }
    return foundcycle;
}

bool shrink2cycles(Network& net, bool unroot) {
    bool foundcycle = false;
    int nh = static_cast<int>(net.hybrid.size());
    int ih = nh;
    while (ih > 0) {
        NodeId h = net.hybrid[static_cast<size_t>(ih - 1)];
        EdgeId minor = net.getparentedgeminor(h);
        EdgeId major = net.getparentedge(h);
        if (net.getparent(minor) != net.getparent(major)) {
            ih -= 1;
            continue;
        }
        foundcycle = true;
        shrink2cycleat(net, minor, major, unroot);
        nh = static_cast<int>(net.hybrid.size());
        ih = nh;
    }
    return foundcycle;
}

// ------------------------------------------------------------------------
// rooting
// ------------------------------------------------------------------------

void rootatnode_name(Network& net, const std::string& name) {
    int found = -1, count = 0;
    for (size_t i = 0; i < net.node.size(); ++i)
        if (net.N(net.node[i]).name == name) {
            if (count == 0) found = static_cast<int>(i);
            ++count;
        }
    if (count == 0) throw SnaqError("node named " + name + " was not found in the network.");
    if (count > 1) throw SnaqError("several nodes were found with name " + name + ".");
    rootatnode_index(net, found);
}

void rootatnode_number(Network& net, int number) {
    int ind = -1;
    for (size_t i = 0; i < net.node.size(); ++i)
        if (net.N(net.node[i]).number == number) {
            ind = static_cast<int>(i);
            break;
        }
    if (ind < 0)
        throw SnaqError("cannot set node " + std::to_string(number) + " as root because it is not part of net");
    rootatnode_index(net, ind);
}

void rootatnode_index(Network& net, int ind) {
    if (ind >= static_cast<int>(net.node.size())) throw SnaqError("node index too large");
    NodeId nd = net.node[static_cast<size_t>(ind)];
    if (net.N(nd).leaf) {
        if (net.N(nd).edge.size() != 1) throw SnaqError("leaf has " + std::to_string(net.N(nd).edge.size()) + " edges!");
        NodeId pn = net.getOtherNode(net.N(nd).edge[0], nd);
        if (net.N(pn).edge.size() <= 2)
            rootatnode_number(net, net.N(pn).number);
        else
            rootonedge_index(net, net.index_of_edge(net.N(nd).edge[0]));
        return;
    }
    int rootsaved = net.rooti;
    net.rooti = ind;
    try {
        directedges(net);
    } catch (const SnaqError& e) {
        // QUIRK: PhyloNetworks converts any error here into a RootMismatch,
        // reverting the root only if the error was a RootMismatch.
        if (dynamic_cast<const RootMismatch*>(&e)) {
            net.rooti = rootsaved;
            directedges(net);
        }
        throw RootMismatch("the desired root is below a reticulation,\nreverting to old root position.");
    }
    if (net.rooti != rootsaved && net.N(net.node[static_cast<size_t>(rootsaved)]).edge.size() == 2)
        fuseedgesat(net, rootsaved);
}

void rootonedge_index(Network& net, int ind) {
    if (ind < 0 || ind >= static_cast<int>(net.edge.size())) throw SnaqError("edge index too large");
    int rootsaved = net.rooti;
    breakedge(net, net.edge[static_cast<size_t>(ind)]);
    net.rooti = static_cast<int>(net.node.size()) - 1;
    try {
        directedges(net);
    } catch (const SnaqError& e) {
        if (dynamic_cast<const RootMismatch*>(&e)) {
            fuseedgesat(net, net.rooti);
            net.rooti = rootsaved;
            directedges(net);
        }
        throw RootMismatch("the desired root is below a reticulation,\nreverting to old root position.");
    }
    if (net.rooti != rootsaved && net.N(net.node[static_cast<size_t>(rootsaved)]).edge.size() == 2)
        fuseedgesat(net, rootsaved);
}

// ------------------------------------------------------------------------
// fliphybrid! (moves_semidirected.jl)
// ------------------------------------------------------------------------

bool fliphybrid(Network& net, NodeId hybridnode, bool minor, bool nohybridladder) {
    bool runDirectEdges = false;
    EdgeId edgetoflip, edgetokeep;
    if (minor) {
        edgetoflip = net.getparentedgeminor(hybridnode);
        edgetokeep = net.getparentedge(hybridnode);
    } else {
        edgetoflip = net.getparentedge(hybridnode);
        edgetokeep = net.getparentedgeminor(hybridnode);
    }
    net.getchildedge(hybridnode);  // errors (as in Julia) if not exactly one child
    NodeId newhybridnode = net.getparent(edgetoflip);
    if (net.N(newhybridnode).hybrid) return false;
    NodeId p2 = net.getparent(edgetokeep);
    std::vector<char> isdesc;
    for (EdgeId e : net.N(newhybridnode).edge) {
        bool isp2desc = false;
        if (e != edgetoflip) {
            NodeId neibr = net.getOtherNode(e, newhybridnode);
            if (!net.N(neibr).leaf && (p2 == neibr || isdescendant_undirected(net, p2, neibr, e))) isp2desc = true;
        }
        isdesc.push_back(isp2desc ? 1 : 0);
    }
    int sum_isdesc = 0;
    for (char c : isdesc) sum_isdesc += c;
    if (sum_isdesc > 1) return false;
    EdgeId newhybridedge;
    if (sum_isdesc == 0) {
        newhybridedge = net.getparentedge(newhybridnode);
    } else {
        size_t k = 0;
        while (!isdesc[k]) ++k;
        newhybridedge = net.N(newhybridnode).edge[k];
        if (net.E(newhybridedge).hybrid) return false;
    }
    if (nohybridladder) {
        if (net.N(net.getOtherNode(newhybridedge, newhybridnode)).hybrid) return false;
        for (EdgeId e : net.N(newhybridnode).edge)
            if (e != edgetoflip && net.E(e).hybrid) return false;
    }
    net.N(hybridnode).hybrid = false;
    net.E(edgetokeep).hybrid = false;
    net.N(newhybridnode).hybrid = true;
    net.E(newhybridedge).hybrid = true;
    net.E(newhybridedge).ismajor = minor;
    net.E(edgetokeep).ismajor = true;
    net.E(edgetoflip).ischild1 = !net.E(edgetoflip).ischild1;
    if (net.getchild(newhybridedge) != newhybridnode) {
        net.E(newhybridedge).ischild1 = !net.E(newhybridedge).ischild1;
        net.rooti = net.index_of_node(hybridnode);
        runDirectEdges = true;
    }
    net.N(newhybridnode).name = net.N(hybridnode).name;
    net.N(hybridnode).name = "";
    net.E(newhybridedge).gamma = net.E(edgetokeep).gamma;
    net.E(edgetokeep).gamma = 1.0;
    auto it = std::find(net.hybrid.begin(), net.hybrid.end(), hybridnode);
    *it = newhybridnode;
    if (runDirectEdges) {
        directedges(net);
    } else {
        for (size_t k = 0; k < net.N(newhybridnode).edge.size(); ++k) {
            EdgeId ce = net.N(newhybridnode).edge[k];
            if (ce == newhybridedge) continue;
            if (net.getparent(ce) != newhybridnode) continue;
            norootbelow(net, ce);
        }
        if (net.E(edgetokeep).containroot) allowrootbelow(net, hybridnode, edgetokeep);
    }
    return true;
}

// ------------------------------------------------------------------------
// addhybridedge! (addHybrid.jl)
// ------------------------------------------------------------------------

std::pair<NodeId, EdgeId> addhybridedge(Network& net, EdgeId edge1, EdgeId edge2, bool hybridpartnernew,
                                        double edgelength, double gamma) {
    if (!(gamma == -1.0 || (gamma <= 1.0 && gamma >= 0.0))) throw SnaqError("invalid γ to add a hybrid edge");
    double gbar = (gamma == -1.0 ? -1.0 : 1.0 - gamma);
    auto [newnode1_tree, edgeabovee1] = breakedge(net, edge1);
    (void)edgeabovee1;
    auto [newnode2_hybrid, edgeabovee2] = breakedge(net, edge2);
    net.N(newnode2_hybrid).hybrid = true;
    net.pushHybrid(newnode2_hybrid);
    EdgeId hybrid_edge = net.new_edge(net.max_edge_number() + 1, edgelength, true, gamma, gamma > 0.5);
    if (hybridpartnernew) {
        net.E(edgeabovee2).hybrid = true;
        net.E(edgeabovee2).gamma = gbar;
        if (gamma > 0.5) net.E(edgeabovee2).ismajor = false;
    } else {
        NodeId c2 = net.getchild(edge2);
        net.rooti = net.index_of_node(c2);
        net.E(edge2).hybrid = true;
        net.E(edge2).gamma = gbar;
        if (gamma > 0.5) net.E(edge2).ismajor = false;
        net.E(edge2).ischild1 = !net.E(edge2).ischild1;
        net.E(edgeabovee2).ischild1 = !net.E(edgeabovee2).ischild1;
    }
    static const std::regex rx("^H(\\d+)$");
    int hnum = net.numhybrids;
    for (NodeId n : net.node) {
        std::smatch m;
        const std::string& nm = net.N(n).name;
        if (std::regex_match(nm, m, rx)) {
            int hi = std::stoi(m[1].str());
            if (!(hnum > hi)) hnum = hi + 1;
        }
    }
    net.N(newnode2_hybrid).name = "H" + std::to_string(hnum);
    net.setNode(hybrid_edge, newnode2_hybrid, newnode1_tree);
    net.setEdge(newnode1_tree, hybrid_edge);
    net.setEdge(newnode2_hybrid, hybrid_edge);
    net.pushEdge(hybrid_edge);
    if (hybridpartnernew) {
        norootbelow(net, edge2);
    } else {
        directedges(net);
        norootbelow(net, edgeabovee2);
    }
    return {newnode2_hybrid, hybrid_edge};
}

bool directionalconflict(const Network& net, NodeId parent, EdgeId edge2, bool hybridpartnernew) {
    if (hybridpartnernew) {
        NodeId c2 = net.getchild(edge2);
        return parent == c2 || isdescendant(net, parent, c2);
    }
    if (!net.E(edge2).containroot || net.E(edge2).hybrid) return true;
    NodeId p2 = net.getparent(edge2);
    return parent == p2 || isdescendant_undirected(net, parent, p2, edge2);
}

// ------------------------------------------------------------------------
// resolvetreepolytomy! (readwrite.jl)
// ------------------------------------------------------------------------

static void resolvetreepolytomy_once(Network& net, NodeId n, EdgeId ce1, EdgeId ce2, bool b1, bool b2) {
    net.removeEdge(n, ce1);
    net.removeEdge(n, ce2);
    net.removeNode(n, ce1);
    net.removeNode(n, ce2);
    // unused_enum = setdiff(1:net.numedges, edge numbers)
    std::set<int> used;
    for (EdgeId e : net.edge) used.insert(net.E(e).number);
    int new_enum = net.numedges + 1;
    for (int k = 1; k <= net.numedges; ++k)
        if (!used.count(k)) {
            new_enum = k;
            break;
        }
    EdgeId e_new = net.new_edge(new_enum, 0.0);
    int max_node = net.max_node_number();
    NodeId n_new = net.new_node(max_node + 1, false, false);
    net.N(n_new).edge.push_back(ce1);
    net.N(n_new).edge.push_back(ce2);
    net.N(n_new).edge.push_back(e_new);
    net.N(n_new).booln1 = net.E(ce1).hybrid || net.E(ce2).hybrid;
    net.setEdge(n, e_new);
    net.setNode(ce1, n_new);
    net.E(ce1).ischild1 = b1;
    net.setNode(ce2, n_new);
    net.E(ce2).ischild1 = b2;
    net.setNode(e_new, n_new, n);
    net.pushNode(n_new);
    net.pushEdge(e_new);
}

void resolvetreepolytomy(Network& net, NodeId n) {
    if (net.N(n).hybrid) throw SnaqError("cannot solve polytomy at a hybrid node " + std::to_string(net.N(n).number) + ".");
    while (net.N(n).edge.size() > 3) {
        // find_twochildedges
        const auto& ee = net.N(n).edge;
        std::vector<char> ischildofn;
        for (EdgeId e : ee) ischildofn.push_back(net.getparent(e) == n);
        int i2 = -1;
        for (int k = static_cast<int>(ee.size()) - 1; k >= 0; --k)
            if (ischildofn[static_cast<size_t>(k)]) {
                i2 = k;
                break;
            }
        EdgeId ce1, ce2;
        bool b1, b2;
        if (i2 < 0) {
            ce1 = ee[0];
            ce2 = ee[1];
            b1 = false;
            b2 = false;
        } else {
            int i1 = -1;
            for (int k = i2 - 1; k >= 0; --k)
                if (ischildofn[static_cast<size_t>(k)]) {
                    i1 = k;
                    break;
                }
            bool notfound = (i1 < 0);
            if (notfound) i1 = (i2 == 0 ? 1 : 0);
            ce1 = ee[static_cast<size_t>(i1)];
            ce2 = ee[static_cast<size_t>(i2)];
            b1 = !notfound;
            b2 = true;
        }
        resolvetreepolytomy_once(net, n, ce1, ce2, b1, b2);
    }
}

// ------------------------------------------------------------------------
// biconnected components (graph_components.jl)
// ------------------------------------------------------------------------

static void bcc_rec(Network& net, NodeId node, int& index, std::vector<EdgeId>& S,
                    std::vector<std::vector<EdgeId>>& blobs, bool ignoreTrivial) {
    net.N(node).intn2 = index;
    net.N(node).intn1 = index;
    ++index;
    for (size_t k = 0; k < net.N(node).edge.size(); ++k) {
        EdgeId e = net.N(node).edge[k];
        NodeId w = (net.E(e).nd[0] == node) ? net.E(e).nd[1] : net.E(e).nd[0];
        if (net.N(w).intn2 == -1) {
            net.E(e).boole2 = true;
            S.push_back(e);
            bcc_rec(net, w, index, S, blobs, ignoreTrivial);
            net.N(node).intn1 = std::min(net.N(node).intn1, net.N(w).intn1);
            if (net.N(w).intn1 >= net.N(node).intn2) {
                std::vector<EdgeId> bb;
                while (!S.empty()) {
                    EdgeId e2 = S.back();
                    S.pop_back();
                    bb.push_back(e2);
                    if (e2 == e) break;
                }
                if (!ignoreTrivial || bb.size() > 1) blobs.push_back(std::move(bb));
            }
        } else if (!net.E(e).boole2 && net.N(node).intn2 > net.N(w).intn2) {
            net.E(e).boole2 = true;
            net.N(node).intn1 = std::min(net.N(node).intn1, net.N(w).intn2);
            S.push_back(e);
        }
    }
}

std::vector<std::vector<EdgeId>> biconnectedcomponents(Network& net, bool ignoreTrivial) {
    for (NodeId n : net.node) {
        net.N(n).intn1 = -1;
        net.N(n).intn2 = -1;
    }
    for (EdgeId e : net.edge) net.E(e).boole2 = false;
    std::vector<EdgeId> S;
    std::vector<std::vector<EdgeId>> blobs;
    int index = 0;
    bcc_rec(net, net.getroot(), index, S, blobs, ignoreTrivial);
    if (!S.empty()) warn("stack of edges not empty at the end of biconnectedcomponents");
    return blobs;
}

std::vector<NodeId> biconnectedcomponent_entrynodes(Network& net, const std::vector<std::vector<EdgeId>>& bcc,
                                                    bool dopreorder) {
    if (dopreorder) {
        directedges(net);
        preorder(net);
    }
    std::vector<NodeId> entry;
    for (const auto& bicomp : bcc) {
        size_t jmin = net.node.size() - 1;  // Julia: jmin = length(net.node) (1-based)
        for (EdgeId e : bicomp) {
            NodeId n = net.getparent(e);
            auto it = std::find(net.vec_node.begin(), net.vec_node.end(), n);
            if (it == net.vec_node.end()) throw SnaqError("node not found in net's pre-ordering 'vec_node'");
            jmin = std::min(jmin, static_cast<size_t>(it - net.vec_node.begin()));
        }
        entry.push_back(net.vec_node[jmin]);
    }
    return entry;
}

std::vector<std::vector<NodeId>> biconnectedcomponent_exitnodes(Network& net,
                                                                const std::vector<std::vector<EdgeId>>& bcc,
                                                                bool dopreorder) {
    if (dopreorder) {
        directedges(net);
        preorder(net);
    }
    std::vector<std::vector<NodeId>> exitnode;
    for (EdgeId e : net.edge) net.E(e).inte1 = -1;
    for (size_t i = 0; i < bcc.size(); ++i)
        for (EdgeId e : bcc[i]) net.E(e).inte1 = static_cast<int>(i) + 1;
    for (size_t i = 0; i < bcc.size(); ++i) {
        std::vector<NodeId> ex;
        for (EdgeId edge : bcc[i]) {
            if (!net.E(edge).ismajor) continue;
            NodeId n = net.getchild(edge);
            for (EdgeId e : net.N(n).edge) {
                if (e == edge) continue;
                if (net.E(e).inte1 != static_cast<int>(i) + 1) {
                    ex.push_back(n);
                    break;
                }
            }
        }
        exitnode.push_back(std::move(ex));
    }
    return exitnode;
}

std::tuple<bool, bool, bool> istreechild(const Network& net) {
    if (net.N(net.getroot()).leaf) throw SnaqError("istreechild requires a root that is not a leaf.");
    std::set<int> hybparent_visited;
    bool hasWatroot = false;
    bool isTCrooted = true;
    for (NodeId h : net.hybrid) {
        for (NodeId n : net.getparents(h)) {
            if (hybparent_visited.count(net.N(n).number)) continue;
            hybparent_visited.insert(net.N(n).number);
            int ntreechild = 0;
            for (EdgeId e : net.N(n).edge)
                if (!net.E(e).hybrid && net.getparent(e) == n) ++ntreechild;
            if (ntreechild > 1 || (net.N(n).hybrid && ntreechild > 0)) continue;
            if (net.N(n).hybrid) return {false, false, false};
            bool atroot = (net.getroot() == n);
            if (ntreechild == 1 && !atroot) continue;
            if (ntreechild == 0) {
                isTCrooted = false;
                if (atroot || !net.E(net.getparentedge(n)).containroot) return {false, false, false};
            }
            if (hasWatroot) {
                if (isTCrooted) throw SnaqError("the network should have already been flagged as not tree-child");
                return {false, false, false};
            }
            hasWatroot = true;
        }
    }
    return {isTCrooted, true, !hasWatroot};
}

}  // namespace csnaq
