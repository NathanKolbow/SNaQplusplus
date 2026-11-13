// Composite log-pseudolikelihood (SNaQ score) and its gradient
// (SNaQ's likelihood/*).
#pragma once

#include <array>
#include <vector>

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
CFMatrix computeexpectedCFmatrix(Network& net, double rho = 0.0, int nthreads = 1);

// computeSNaQscore!(N, q, ρ): score of a copy of N made semidirected, on all
// quartets. Also stored in N.fscore.
double computeSNaQscore_net(Network& N, const CFMatrix& q, double rho = 0.0, int nthreads = 1);

// computeexpectedCF(qdata, params, ρ) -> (eCF1, eCF2)
std::array<double, 2> computeexpectedCF(const QuartetDataImpl& qd, const std::vector<double>& params, double rho);

}  // namespace csnaq
