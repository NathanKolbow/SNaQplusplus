// Phylogenetic network data structure mirroring PhyloNetworks.jl's
// `HybridNetwork`, `Node` and `Edge`.
//
// Nodes and edges live in pools and are referred to by integer ids. The
// network's ordered lists (`node`, `edge`, `leaf`, `hybrid`) hold ids in the
// same order as PhyloNetworks' vectors: SNaQ's behaviour depends on these
// orders (parameter numbering, traversal order, which edge is fused, ...),
// so every operation below reproduces PhyloNetworks' list manipulations.
#pragma once

#include <string>
#include <vector>

#include "util.hpp"

namespace csnaq {

using NodeId = int;
using EdgeId = int;
constexpr int NONE = -1;

// Small list of nodes (children / parents): no heap allocation for degree <= 4.
using NodeList = SmallVec<int, 4>;

struct Node {
    int number = -1;
    bool leaf = false;
    bool hybrid = false;
    std::string name;
    SmallVec<EdgeId, 4> edge;  // incident edges, in PhyloNetworks' order
    // scratch attributes used by some algorithms
    bool booln1 = false;  // incident to a hybrid edge?
    int intn1 = -1;       // lowpoint (biconnected components)
    int intn2 = -1;       // discovery index (biconnected components)
};

struct Edge {
    int number = -1;
    double length = 1.0;
    bool hybrid = false;
    double y = -1.0;
    double z = -1.0;
    double gamma = 1.0;
    NodeId nd[2] = {NONE, NONE};
    int nnode = 0;  // number of attached nodes (0..2)
    bool ischild1 = true;
    bool ismajor = true;
    int inte1 = -1;
    bool containroot = true;
    bool boole1 = true;
    bool boole2 = false;

    NodeId node(int i) const { return nd[i]; }  // 0-based (Julia's node[i+1])
    void push_node(NodeId n) {
        if (nnode >= 2) throw SnaqError("vector of nodes already has 2 values");
        nd[nnode++] = n;
    }
    void remove_node(NodeId n);  // Julia removeNode!(n, edge)
    bool has_node(NodeId n) const { return (nnode > 0 && nd[0] == n) || (nnode > 1 && nd[1] == n); }
};

class Network {
public:
    // pools
    std::vector<Node> npool;
    std::vector<Edge> epool;

    // PhyloNetworks lists (ids into the pools)
    std::vector<NodeId> node;
    std::vector<EdgeId> edge;
    std::vector<NodeId> leaf;
    std::vector<NodeId> hybrid;

    int numtaxa = 0;
    int numnodes = 0;
    int numedges = 0;
    int numhybrids = 0;
    int rooti = -1;  // 0-based index into `node`
    std::vector<std::string> names;
    double fscore = 0.0;
    bool isrooted = false;
    std::vector<NodeId> vec_node;  // preorder

    // ---- pool access ----
    Node& N(NodeId id) { return npool[static_cast<size_t>(id)]; }
    const Node& N(NodeId id) const { return npool[static_cast<size_t>(id)]; }
    Edge& E(EdgeId id) { return epool[static_cast<size_t>(id)]; }
    const Edge& E(EdgeId id) const { return epool[static_cast<size_t>(id)]; }

    NodeId new_node(int number, bool leaf, bool hybrid = false);
    EdgeId new_edge(int number, double length = 1.0);
    // Julia: Edge(number, length, hybrid, gamma, ismajor)
    EdgeId new_edge(int number, double length, bool hybrid, double gamma, bool ismajor);

    // ---- basic queries (PhyloNetworks auxiliary.jl) ----
    NodeId getchild(EdgeId e) const { return E(e).nd[E(e).ischild1 ? 0 : 1]; }
    NodeId getparent(EdgeId e) const { return E(e).nd[E(e).ischild1 ? 1 : 0]; }
    NodeId getOtherNode(EdgeId e, NodeId n) const { return E(e).nd[0] == n ? E(e).nd[1] : E(e).nd[0]; }
    bool isparentof(NodeId n, EdgeId e) const { return n == getparent(e); }
    bool ischildof(NodeId n, EdgeId e) const { return n == getchild(e); }
    NodeList getchildren(NodeId n) const;
    NodeList getparents(NodeId n) const;
    NodeId getparent_node(NodeId n) const;  // major parent node
    EdgeId getparentedge(NodeId n) const;   // major parent edge
    EdgeId getparentedgeminor(NodeId n) const;
    EdgeId getchildedge(NodeId n) const;
    NodeId getchild_node(NodeId n) const { return getchild(getchildedge(n)); }
    EdgeId getpartneredge(EdgeId e, NodeId n) const;
    EdgeId getpartneredge(EdgeId e) const { return getpartneredge(e, getchild(e)); }
    EdgeId getconnectingedge(NodeId a, NodeId b) const;  // NONE if not connected
    NodeId getroot() const { return node[static_cast<size_t>(rooti)]; }
    bool isrootof(NodeId n) const { return n == getroot(); }
    bool isexternal(EdgeId e) const;
    bool isTree() const;
    int index_of_node(NodeId n) const;  // index in `node`, -1 if absent
    int index_of_edge(EdgeId e) const;  // index in `edge`, -1 if absent
    NodeId find_node_by_number(int number) const;  // NONE if absent
    std::vector<std::string> tiplabels() const;
    int max_edge_number() const;
    int max_node_number() const;

    // ---- list maintenance (PhyloNetworks auxiliary.jl) ----
    void pushEdge(EdgeId e);
    void pushNode(NodeId n);
    void pushHybrid(NodeId n);
    void setNode(EdgeId e, NodeId n);               // setNode!(edge, node)
    void setNode(EdgeId e, NodeId n1, NodeId n2);   // setNode!(edge, [n1,n2])
    void setEdge(NodeId n, EdgeId e);               // setEdge!(node, edge)
    void removeEdge(NodeId n, EdgeId e);            // removeEdge!(node, edge)
    void removeNode(NodeId n, EdgeId e);            // removeNode!(node, edge)
    void deleteEdge(EdgeId e);                      // deleteEdge!(net, e; part=false)
    void deleteNode(NodeId n);                      // deleteNode!(net, n)
    void removeHybrid(NodeId n);                    // removeHybrid!(net, n)
    void removeLeaf(NodeId n);                      // removeLeaf!(net, n)

    // ---- copies ----
    // Base.deepcopy: exact copy (pools are value types).
    Network clone() const { return *this; }
    // SNaQ's `deepcopynetwork`: compacts pools, rebuilds each node's edge list
    // in `net.edge` order, and resets attributes it does not copy.
    Network deepcopynetwork() const;
};

}  // namespace csnaq
