// Expected-CF equations of quartet classes whose four lineages meet through
// reticulations, built without copying or editing the network.
//
// SNaQ.jl (and portreticulateequation, kept as the reference implementation)
// finds the lowest hybrid node of the quartet's 4-taxon network, copies the
// network once per way the lineages below it can pass it, edits each copy,
// and recurses. Here the network never changes: the state of the recursion is
// the position of the four lineages, which only move up.
//
//  - Each lineage climbs through tree nodes until it reaches a hybrid node it
//    has not passed yet, or `top`. Lineages that reach the same end node meet
//    on the way, along tree-like paths. (The paths are not walked: the nodes
//    of a path are the ancestors of its start that are not above its end, and
//    the edges that two lineages share are counted and listed from the
//    intersection of two such sets.)
//  - Four lineages with the same end, or three that the fourth cannot join
//    before they have all met: the equation is that of a tree. (Three
//    lineages that reach a node without having coalesced are exchangeable, so
//    nothing above that node matters.)
//  - Otherwise a lowest end hybrid node H is passed (no lineage is below
//    another end hybrid that is below H):
//      one lineage at H: it takes the minor or the major edge (2 divisions),
//        unless both ways lead it, alone, to the same node: it is then moved
//        there and H does not appear in the equation;
//      two lineages at H: they may coalesce on the path they share below H,
//        then each takes the minor or the major edge (4 divisions).
//
// The result has the same value as the reference for every parameter value.
// It can be a smaller tree: the reference removes reticulations that cannot
// matter only in the two simplest cases.
#pragma once

#include <array>
#include <cstdint>
#include <vector>

#include "equations.hpp"
#include "network.hpp"

namespace snaqpp {

class PortEquationBuilder {
public:
    // `net` must be directed consistently from its root, and must not change
    // while the builder is used. Parameters of the equations are pids.
    explicit PortEquationBuilder(const Network& net);

    // Equation of four lineages that enter the network at the nodes `ports`,
    // in the order (1, 2, 3, 4): which_coal = 1 means 12|34, 2 means 13|24,
    // 3 means 14|23. `top` is a node that every path from the root to a port
    // goes through. Thread-safe.
    RecursiveCFEquation build(const std::array<NodeId, 4>& ports, NodeId top) const;

private:
    // A lineage is at `node`, about to climb. `via`: the parent edge of the
    // hybrid node `node` that it takes (NONE: not chosen yet, or tree node).
    struct Pos {
        NodeId node;
        EdgeId via;
    };
    using State = std::array<Pos, 4>;
    RecursiveCFEquation recur(const State& st, NodeId top, int depth) const;
    const uint64_t* ancestors(const Pos& p) const {
        const NodeId n = p.via == NONE ? p.node : eparent_[static_cast<size_t>(p.via)];
        return &anc_[static_cast<size_t>(npos_[static_cast<size_t>(n)]) * W_];
    }
    bool isancestor(NodeId a, const uint64_t* set) const {
        const size_t b = static_cast<size_t>(npos_[static_cast<size_t>(a)]);
        return (set[b >> 6] >> (b & 63)) & 1;
    }

    // by NodeId
    std::vector<EdgeId> pmaj_, pmin_;  // parent edges; pmin_ != NONE: hybrid node
    std::vector<int> npid_, npos_;     // pid; position in net.node
    // where a lineage at the node stops climbing through tree nodes: the first
    // hybrid node at or above it, or the root
    std::vector<NodeId> tend_;
    std::vector<NodeId> bypos_;  // by position in net.node: the node
    // by EdgeId
    std::vector<NodeId> eparent_, echild_;
    std::vector<int> epid_;
    // by position in net.node: the node and its ancestors (bitset over positions)
    size_t W_ = 0;
    std::vector<uint64_t> anc_;
    int maxdepth_ = 0;
};

// Number of nodes of an equation tree.
size_t equationsize(const RecursiveCFEquation& eqn);

}  // namespace snaqpp
