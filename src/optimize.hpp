// Numerical optimization of branch lengths and γ (SNaQ's fit_parameters.jl):
// with NLopt's LD_LBFGS exactly as SNaQ.jl does, or with a projected Newton
// method that uses the curvature of the quartet-class likelihood.
#pragma once

#include <vector>

#include "equations.hpp"
#include "likelihood.hpp"

namespace csnaq {

enum class FitMethod {
    lbfgs,   // NLopt's LD_LBFGS, as SNaQ.jl
    newton,  // projected Newton with the curvature of evalclasses (quartet classes only)
};

struct FitOptions {
    // Maximum number of evaluations of the score.
    // -1: max(30, numtaxa) with lbfgs, 200 with newton.
    int maxeval = -1;
    double ftolRel = 1e-12;
    double ftolAbs = 1e-12;
    double xtolRel = 1e-12;
    double xtolAbs = 1e-12;
    FitMethod method = FitMethod::lbfgs;
    // Shallow evaluation of reticulate classes during the fit (fitclasses
    // only; 0: exact). Sub-equations whose weight at the starting parameters
    // is below `shallow` are approximated as `shallowmode` says (see
    // EqnDag::shallowplan). The score returned is the exact score at the
    // fitted parameters.
    double shallow = 0.0;
    ShallowMode shallowmode = ShallowMode::path;
};

// fitnumericalparameters!(net, eqns, observed_CFs, ρ; maxeval, ftolRel, ...)
// Returns the optimized SNaQ score, also stored in net.fscore.
double fitnumericalparameters(Network& net, const std::vector<QuartetData>& eqns, const CFMatrix& obs, double rho,
                              const FitOptions& opt, int nthreads = 1);

// fitnumericalparameters!(net, oCFs): all quartets, equations built here.
double fitnumericalparameters_all(Network& net, const CFMatrix& obs, double rho, const FitOptions& opt,
                                  int nthreads = 1);

}  // namespace csnaq
