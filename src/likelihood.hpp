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

// Quartets that no reticulation divides, flattened (TreeQuartetBatch).
struct TreeQuartetBatch {
    std::vector<int> coal;
    std::vector<size_t> off;
    std::vector<double> qres, qoth, qlogq;
    std::vector<size_t> general;  // indices of quartets a hybrid divides
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

// computeexpectedCF(qdata, params, ρ) -> (eCF1, eCF2)
std::array<double, 2> computeexpectedCF(const QuartetDataImpl& qd, const std::vector<double>& params, double rho);

}  // namespace csnaq
