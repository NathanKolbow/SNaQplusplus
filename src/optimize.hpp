// Numerical optimization of branch lengths and γ (SNaQ's fit_parameters.jl),
// using NLopt's LD_LBFGS exactly as SNaQ.jl does.
#pragma once

#include <vector>

#include "equations.hpp"
#include "likelihood.hpp"

namespace csnaq {

struct FitOptions {
    int maxeval = -1;  // -1: max(30, numtaxa)
    double ftolRel = 1e-12;
    double ftolAbs = 1e-12;
    double xtolRel = 1e-12;
    double xtolAbs = 1e-12;
};

// fitnumericalparameters!(net, eqns, observed_CFs, ρ; maxeval, ftolRel, ...)
// Returns the optimized SNaQ score, also stored in net.fscore.
double fitnumericalparameters(Network& net, const std::vector<QuartetData>& eqns, const CFMatrix& obs, double rho,
                              const FitOptions& opt, int nthreads = 1);

// fitnumericalparameters!(net, oCFs): all quartets, equations built here.
double fitnumericalparameters_all(Network& net, const CFMatrix& obs, double rho, const FitOptions& opt,
                                  int nthreads = 1);

}  // namespace csnaq
