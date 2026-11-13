// The equations of the reticulate classes of a network, with every distinct
// sub-equation stored once.
//
// The equation of a class is a tree (RecursiveCFEquation). Many of its
// sub-equations appear several times in it, and many more in the equations of
// other classes: once their lineages have passed the reticulations that tell
// them apart, two classes continue the same way. EqnDag stores each distinct
// sub-equation of a set of classes as one node, children before parents, so
// that one evaluation computes the value and the derivatives of each of them
// once, in two loops over flat arrays:
//
//  - values:      exp(-sum of lengths) of each distinct set of edges, then the
//                 expected CFs (e1, e2) of each node from those of its children;
//  - derivatives: d(e1, e2) / d(parameter) of each node, for each of its
//                 parameters, from those of its children.
//
// A class then reads the expected CFs and derivatives of its root node.
#pragma once

#include <cstdint>
#include <vector>

#include "equations.hpp"

namespace csnaq {

struct ClassEqn;

// 128-bit content hash of a sub-equation, and its number of nodes.
struct EqnHash {
    uint64_t a = 0, b = 0;
    uint32_t size = 0;
};

// Hashes of the nodes of `eqn`, in preorder. Sub-equations with the same hash
// have the same value for every parameter value (the order of the edges of a
// node, and of the two divisions of a 4-way division where one lineage takes
// the minor edge and the other the major edge, do not count).
std::vector<EqnHash> hashequation(const RecursiveCFEquation& eqn);

// What replaces a sub-equation of low weight in a shallow evaluation
// (EqnDag::shallowplan).
enum class ShallowMode {
    path,     // a single way: the lineages take the edge of largest γ
    average,  // the average of the tree-like equations it leads to
    freeze,   // its value at the parameters the plan was made at
};

class EqnDag {
public:
    struct Node {
        uint8_t kind = 0;   // 0: tree-like (leaf of the equation); 2, 4: number of divisions
        uint8_t wc = 0;     // which_coal (4 divisions: 0 if the lineages cannot coalesce here)
        int32_t H = -1;     // pid of γ
        int32_t set = -1;   // its set of edges (-1: none)
        uint32_t child[4] = {0, 0, 0, 0};
        // Its parameters. Divisions: pid[poff .. poff + np), sorted. Leaves:
        // the edges of its set (poff is not used).
        uint32_t poff = 0, np = 0;
        uint32_t hpos = 0;  // position of H among its parameters
        // map[moff ..]: the position, among its parameters, of each parameter
        // of child 0, then of child 1, ..., then of each edge of its set
        uint32_t moff = 0;
    };

    // Values of one evaluation.
    struct Work {
        std::vector<double> ex;   // per set: exp(-sum of lengths)
        std::vector<double> val;  // per node: e1, e2
        std::vector<double> jac;  // per parameter of each division node (as in pid): d e1, d e2
    };

    // An approximate ("shallow") evaluation: what to do with each node.
    struct Plan {
        enum : uint8_t { exact = 0, skip = 1, path = 2, average = 3, frozen = 4 };
        std::vector<uint8_t> act;
        std::vector<uint8_t> way;       // path: the division that is followed
        std::vector<double> frozenval;  // frozen: e1, e2
        // average: the leaves under node i are avgleaf[avgoff[i] .. avgoff[i] + avgcnt[i]);
        // avgpos lists, for each of them in turn, the position of each of its
        // edges among the parameters of node i, from avgposoff[i]
        std::vector<uint32_t> avgoff, avgcnt, avgposoff, avgleaf, avgpos;
        size_t nshallow = 0;  // division nodes that are approximated
        size_t nskipped = 0;  // division nodes that are not evaluated at all
        size_t ndivisions = 0;
    };

    // Stores the sub-equations of `classes`; root[c] is the node of classes[c].
    void build(const std::vector<const ClassEqn*>& classes);
    bool empty() const { return node.empty(); }

    // With a plan: the approximate evaluation it describes.
    void values(const double* params, double alpha, Work& w, const Plan* plan = nullptr) const;
    // After values(params, alpha, w, plan).
    void derivatives(const double* params, double alpha, Work& w, const Plan* plan = nullptr) const;

    // The weight of a node is the largest product, over the ways from the root
    // of a class down to it, of the probabilities on the way (γ, 1 - γ, the
    // probability of not coalescing before a hybrid node, ...), at `params`:
    // an upper bound of what the node contributes to any expected CF.
    // The plan approximates the division nodes of weight below `eps` as `mode`
    // says, and skips the nodes that are then no longer needed. Which nodes
    // are approximated does not change with the parameters afterwards, so the
    // approximate score is a smooth function.
    Plan shallowplan(const double* params, double alpha, double eps, ShallowMode mode, Work& w) const;

    std::vector<Node> node;       // children before parents
    std::vector<uint32_t> setoff;  // edges of set s: setpid[setoff[s] .. setoff[s + 1])
    std::vector<int> setpid;
    std::vector<int> pid;
    std::vector<uint32_t> map;
    std::vector<uint32_t> root;
};

}  // namespace csnaq
