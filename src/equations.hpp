// Expected-CF equations for each quartet of a network (SNaQ's
// quartet_equations/*), and the mapping between network objects and the
// numerical parameters that are optimized (optimization/parameters.jl).
#pragma once

#include <array>
#include <memory>
#include <string>
#include <vector>

#include "network.hpp"

namespace csnaq {

// Number -> parameter index map (SNaQ's `param_map::Dict{Int,Int}`), keyed by
// edge numbers (branch lengths) and hybrid node numbers (γ).
class ParamMap {
public:
    void clear() {
        base_ = 0;
        v_.clear();
    }
    void set(int number, int idx);
    bool has(int number) const {
        long long k = static_cast<long long>(number) - base_;
        return k >= 0 && k < static_cast<long long>(v_.size()) && v_[static_cast<size_t>(k)] >= 0;
    }
    int at(int number) const;  // throws if missing (Julia KeyError)
    size_t size() const { return count_; }

private:
    int base_ = 0;
    std::vector<int> v_;
    size_t count_ = 0;
};

// An optimized parameter is either an edge's length or a hybrid node's minor γ.
struct ParamObj {
    bool is_hybrid_node;  // true: γ of hybrid node `id`; false: length of edge `id`
    int id;
};

struct OptInfo {
    int narg = 0;
    ParamMap param_map;
    std::vector<ParamObj> idx_obj;  // parameter index -> object
    std::vector<double> params;
    std::vector<double> LB, UB, init_steps;
};

// gatheroptimizationinfo(net, change_numbers)
OptInfo gatheroptimizationinfo(Network& net, bool change_numbers = true);
// countoptimizationparams(net)
int countoptimizationparams(const Network& net);
// gatherparams(net, param_map)
std::vector<double> gatherparams(const Network& net, const ParamMap& pm);
// setX!(net, X, idx_obj_map)
void setX(Network& net, const double* X, const std::vector<ParamObj>& idx_obj);

struct RecursiveCFEquation {
    bool can_coalesce_here = true;
    std::vector<int> coal_edges;  // 0-based parameter indices
    int which_coal = 0;           // 0 = NA, 1 = ab|cd, 2 = ac|bd, 3 = ad|bc
    int division_H = -1;          // γ parameter index; -1: not set
    std::vector<RecursiveCFEquation> divisions;
};

struct QuartetDataImpl {
    RecursiveCFEquation eqn;
    std::vector<int> relevant_params;
    std::array<int, 4> q_taxa;  // 0-based indices into the sorted taxon list
};
// Shared and immutable: unchanged quartets are reused across networks,
// as SNaQ.jl does when it copies `old_eqns[j]` into the new equations.
using QuartetData = std::shared_ptr<const QuartetDataImpl>;

// Per-network data reused across quartets (SNaQ's TreeQuartetContext).
struct TreeQuartetContext {
    const Network* net = nullptr;
    std::vector<std::string> taxa;           // sorted taxon names
    std::vector<NodeId> leaf_of_taxon;       // taxon index -> leaf node
    std::vector<std::vector<NodeId>> chains; // taxon index -> ancestor chain
    bool hashybrids = false;
    std::vector<int> nodepos;  // NodeId -> position in net.node
    std::vector<int> edgepos;  // EdgeId -> position in net.edge
};
TreeQuartetContext treequartetcontext(const Network& net);

// fixnegativeedges!(net)
void fixnegativeedges(Network& net, bool shouldwarn = true);

struct EquationSet {
    std::vector<QuartetData> eqns;
    OptInfo info;  // from gatheroptimizationinfo(net) (change_numbers = true)
    std::vector<std::string> taxa;
};

// findquartetequations(net, sampled_quartets): `ranks` are 0-based.
EquationSet findquartetequations(Network& net, const std::vector<long long>& ranks, int nthreads = 1);
// All quartets.
EquationSet findquartetequations(Network& net, int nthreads = 1);

// rebuildquartetequations!(net, old_eqns, new_eqns)
void rebuildquartetequations(Network& net, const std::vector<QuartetData>& old_eqns,
                             std::vector<QuartetData>& new_eqns, int nthreads = 1);

// findquartetequations4taxa(ctx, taxa, param_map)
QuartetData findquartetequations4taxa(const TreeQuartetContext& ctx, const std::array<int, 4>& taxa,
                                      const ParamMap& pm);

// In-place updates after rNNI(1) / rNNI(2) moves (inplace_updates.jl).
void applyrNNI1update(Network& Nprime, const std::vector<QuartetData>& old_q, std::vector<QuartetData>& new_q,
                      const ParamMap& pm, NodeId u, int nthreads);
void applyrNNI2update(Network& Nprime, const std::vector<QuartetData>& old_q, std::vector<QuartetData>& new_q,
                      const ParamMap& pm, NodeId s, NodeId t, NodeId u, NodeId v, int nthreads);

}  // namespace csnaq
