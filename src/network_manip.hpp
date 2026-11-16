// Ports of the PhyloNetworks.jl network manipulation functions used by SNaQ.
#pragma once

#include <tuple>
#include <vector>

#include "network.hpp"

namespace snaqpp {

struct DeleteOpts {
    bool nofuse = false;
    bool simplify = true;
    bool unroot = false;
    bool multgammas = false;
    bool keeporiginalroot = false;
};

inline double addBL(double x, double y) { return (x == -1.0 || y == -1.0) ? -1.0 : x + y; }
inline double multiplygammas(double x, double y) { return (x == -1.0 || y == -1.0) ? -1.0 : x * y; }

// directedges!(net; checkMajor)
void directedges(Network& net, bool checkMajor = true);
// preorder!(net): fills net.vec_node
void preorder(Network& net);

void allowrootbelow(Network& net, EdgeId e);
void allowrootbelow(Network& net, NodeId n, EdgeId pe);
void allowrootbelow(Network& net);
void norootbelow(Network& net, EdgeId e);

bool isdescendant(const Network& net, NodeId des, NodeId anc);
bool isdescendant_undirected(const Network& net, NodeId des, NodeId anc, EdgeId parentedge);

// deletehybridedge!(net, edge, nofuse, unroot, multgammas, simplify, keeporiginalroot)
void deletehybridedge(Network& net, EdgeId edge, bool nofuse = false, bool unroot = false, bool multgammas = false,
                      bool simplify = true, bool keeporiginalroot = false);
// deleteleaf!(net, i; index=true, ...): i is a 0-based index into net.node
void deleteleaf_index(Network& net, int i, const DeleteOpts& o);
// deleteleaf!(net, number; index=false, ...)
void deleteleaf_number(Network& net, int nodeNumber, const DeleteOpts& o);
// deleteleaf!(net, node::Node; ...)
void deleteleaf_node(Network& net, NodeId n, const DeleteOpts& o);
// deleteleaf!(net, name::String; ...)
void deleteleaf_name(Network& net, const std::string& name, const DeleteOpts& o);

NodeId addleaf(Network& net, NodeId speciesnode, const std::string& leafname, double edgelength = -1.0);
NodeId addleaf_on_edge(Network& net, EdgeId startingedge, const std::string& leafname, double edgelength = -1.0);

// breakedge!(edge, net) -> (newnode, newedge)
std::pair<NodeId, EdgeId> breakedge(Network& net, EdgeId edge);
// fuseedgesat!(i, net, multgammas): i is a 0-based index into net.node
EdgeId fuseedgesat(Network& net, int i, bool multgammas = false);
void removedegree2nodes(Network& net, bool keeproot = false);
bool shrink2cycles(Network& net, bool unroot = false);
bool shrink3cycles(Network& net, bool unroot = false);

void rootatnode_index(Network& net, int ind);
void rootatnode_number(Network& net, int number);
void rootatnode_name(Network& net, const std::string& name);
void rootonedge_index(Network& net, int ind);

// fliphybrid!(net, hybridnode, minor=true, nohybridladder=false).
// Returns false if the flip was not possible (network unchanged).
bool fliphybrid(Network& net, NodeId hybridnode, bool minor = true, bool nohybridladder = false);

// addhybridedge!(net, edge1, edge2, hybridpartnernew, edgelength, gamma)
std::pair<NodeId, EdgeId> addhybridedge(Network& net, EdgeId edge1, EdgeId edge2, bool hybridpartnernew,
                                        double edgelength = -1.0, double gamma = -1.0);
bool directionalconflict(const Network& net, NodeId parent, EdgeId edge2, bool hybridpartnernew);

void resolvetreepolytomy(Network& net, NodeId n);

// biconnectedcomponents(net, ignoreTrivial)
std::vector<std::vector<EdgeId>> biconnectedcomponents(Network& net, bool ignoreTrivial = false);
std::vector<NodeId> biconnectedcomponent_entrynodes(Network& net, const std::vector<std::vector<EdgeId>>& bcc,
                                                    bool preorder = true);
std::vector<std::vector<NodeId>> biconnectedcomponent_exitnodes(Network& net,
                                                                const std::vector<std::vector<EdgeId>>& bcc,
                                                                bool preorder = true);

// istreechild(net) -> (rooted tree-child, weakly tree-child, strongly tree-child)
std::tuple<bool, bool, bool> istreechild(const Network& net);

}  // namespace snaqpp
