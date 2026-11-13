// Composite log-pseudolikelihood (SNaQ score) and its gradient
// (SNaQ's likelihood/*).
#pragma once

#include <array>
#include <vector>

#include "eqndag.hpp"
#include "equations.hpp"

namespace csnaq {

// Observed CFs: row-major nq x 3.
struct CFMatrix {
    std::vector<double> v;
    size_t nrow() const { return v.size() / 3; }
    double operator()(size_t i, int k) const { return v[3 * i + static_cast<size_t>(k)]; }
    double& operator()(size_t i, int k) { return v[3 * i + static_cast<size_t>(k)]; }
    const double* row(size_t i) const { return &v[3 * i]; }
};

double rhotoalpha(double rho);
double alphatorho(double alpha);

// The quartets of a fit, prepared for repeated loss/gradient evaluation
// (SNaQ's TreeQuartetBatch, extended).
//
// Quartets are merged into groups that have identical expected CFs as a
// function of the parameters: tree-like quartets that coalesce over the same
// set of branches, and reticulate quartets with identical equation trees.
// A group's contribution to the loss is linear in its quartets' observed CFs,
// so a group is evaluated once with the sum of their observed CFs. The result
// equals SNaQ.jl's per-quartet sum up to floating-point summation order, and
// the work per evaluation drops from O(#quartets) to O(#groups).
struct TreeQuartetBatch {
    // tree-like groups
    std::vector<int> coal;
    std::vector<size_t> off;
    std::vector<double> qres, qoth, qlogq;
    // reticulate groups
    std::vector<const QuartetDataImpl*> gen_rep;  // representative equation
    std::vector<std::array<double, 3>> gen_q;     // summed observed CFs
    std::vector<double> gen_qlogq;                // sum of q log q over the group
    std::vector<std::vector<int>> gen_params;     // parameters in the equation tree
    size_t nquartets = 0;
};
TreeQuartetBatch treequartetbatch(const std::vector<QuartetData>& qdata, const CFMatrix& obs);

// Per-thread scratch for the general (reticulate) quartets.
class LossWorkspace {
public:
    LossWorkspace() = default;
    struct Impl;
    ~LossWorkspace();
    LossWorkspace(const LossWorkspace&) = delete;
    LossWorkspace& operator=(const LossWorkspace&) = delete;
    Impl* get(int nparam, int nthreads);

private:
    Impl* impl_ = nullptr;
};

// ------------------------------------------------------------------------
// Quartet classes
// ------------------------------------------------------------------------
// A class is the set of quartets whose four lineages enter the network
// through the same four cut edges: they have the same expected CFs, so the
// fit needs one equation per class and the sum of the class's observed CFs.
// Class equations use persistent parameter ids (Node::pid / Edge::pid), so
// they stay valid in copies of the network and after topology moves.

struct ClassEqn {
    RecursiveCFEquation eqn;    // parameters are pids
    std::vector<int> params;    // the distinct pids of the equation tree
    std::vector<EqnHash> hash;  // hashequation(eqn), for reticulate classes
};

// The classes of a network, prepared for repeated loss/gradient evaluation.
struct ClassBatch {
    int npid = 0;  // parameter arrays are indexed by pid: 0 .. npid-1
    // tree-like classes (off.size() = number of classes + 1)
    std::vector<int> coal;
    std::vector<uint32_t> off;
    std::vector<double> qres, qoth;
    // reticulate classes
    std::vector<const ClassEqn*> gen;
    std::vector<std::array<double, 3>> genQ;  // summed observed CFs, in the class's column order
    // The equations of `gen` with their sub-equations shared (see preparebatch).
    // Classes with the same equation are evaluated together: `root` lists the
    // distinct root nodes, `rootQ` the summed observed CFs of their classes.
    EqnDag dag;
    std::vector<uint32_t> root;
    std::vector<std::array<double, 3>> rootQ;
    double qlogq = 0.0;  // sum over all quartets of q log q
};
// Fills batch.dag, root and rootQ from batch.gen and genQ.
void preparebatch(ClassBatch& batch);

// Loss and its gradient with respect to every pid (`params` and `grad` have
// `npid` entries). If `hess` is given, it receives a positive semidefinite
// approximation of minus the Hessian over the parameters with dense[pid] >= 0
// (ndense x ndense, row-major): exact for tree-like classes, Gauss-Newton for
// reticulate classes.
double evalclasses(const ClassBatch& batch, const double* params, double* grad, double alpha, LossWorkspace& ws,
                   int nthreads, const int* dense = nullptr, int ndense = 0, double* hess = nullptr);
// The two halves of evalclasses. evalclassvalues returns the loss and keeps
// what the derivatives need in `ws`; evalclassderivatives must follow a call
// of evalclassvalues with the same batch, params, alpha and ws.
// With `plan` (EqnDag::shallowplan of batch.dag; the same in both calls): the
// approximate evaluation it describes.
double evalclassvalues(const ClassBatch& batch, const double* params, double alpha, LossWorkspace& ws,
                       int nthreads, const EqnDag::Plan* plan = nullptr);
void evalclassderivatives(const ClassBatch& batch, const double* params, double* grad, double alpha,
                          LossWorkspace& ws, int nthreads, const int* dense = nullptr, int ndense = 0,
                          double* hess = nullptr, const EqnDag::Plan* plan = nullptr);
// evalclasses with each reticulate class evaluated on its own, by recursion on
// its equation tree, and the curvature summed term by term: for validation.
double evalclasses_reference(const ClassBatch& batch, const double* params, double* grad, double alpha,
                             LossWorkspace& ws, int nthreads, const int* dense = nullptr, int ndense = 0,
                             double* hess = nullptr);

// Probability that two lineages at a hybrid node with minor γ `g` both take
// the minor edge (type 1), both the major edge (2), or one each (3, 4), with
// inheritance correlation 1 / (1 + a); and its derivative in g.
double quadsplitprobability(int type, double g, double a);
double quadsplitprobabilitygradient(int type, double g, double a);

// Expected CFs (in the class's column order) of one class.
std::array<double, 3> classexpectedCF(const ClassEqn& ce, const double* params, int npid, double alpha);

// Upper bounds of one class's term of the loss, over all parameter values
// (the class is given the values that suit it best, whatever the other
// classes need). Tree-like class: `qres` is the summed observed CF of the
// split the network displays, `qoth` of the other two, and `hasparams` tells
// whether its path carries any parameter. General class: summed observed CFs.
double treeclassbound(double qres, double qoth, bool hasparams);
double generalclassbound(const std::array<double, 3>& Q);

// computelossandgradient!(batch, qdata, params, grad, q, α)
double computelossandgradient(const TreeQuartetBatch& batch, const std::vector<QuartetData>& qdata,
                              const double* params, size_t nparam, double* grad, const CFMatrix& q, double alpha,
                              LossWorkspace& ws, int nthreads);
// computelossandgradient!(qdata, params, grad, q, α) (no batch)
double computelossandgradient(const std::vector<QuartetData>& qdata, const std::vector<double>& params,
                              std::vector<double>& grad, const CFMatrix& q, double alpha, int nthreads = 1);

// computeSNaQscore!(qdata, params, q, ρ)
double computeSNaQscore(const std::vector<QuartetData>& qdata, const std::vector<double>& params, const CFMatrix& q,
                        double rho, int nthreads = 1);

// computeexpectedCFmatrix(net, ρ): expected CFs of all quartets of `net`
// (modifies net as SNaQ.jl does: renumbering, negative hybrid lengths -> 0).
// Computed from the network's quartet classes.
CFMatrix computeexpectedCFmatrix(Network& net, double rho = 0.0, int nthreads = 1);

// computeSNaQscore!(N, q, ρ): score of a copy of N made semidirected, on all
// quartets. Also stored in N.fscore.
double computeSNaQscore_net(Network& N, const CFMatrix& q, double rho = 0.0, int nthreads = 1);

// computeexpectedCF(qdata, params, ρ) -> (eCF1, eCF2)
std::array<double, 2> computeexpectedCF(const QuartetDataImpl& qd, const std::vector<double>& params, double rho);

}  // namespace csnaq
