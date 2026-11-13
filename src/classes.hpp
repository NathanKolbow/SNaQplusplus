// Quartet classes: the sufficient statistics of a network's pseudolikelihood.
//
// Two quartets are in the same class when their four lineages enter the
// network through the same four cut edges ("ports": for each taxon, the
// highest cut edge above it whose descendants include none of the other three
// taxa). Below its port a lineage is alone, so nothing there affects the
// quartet's expected CFs: all quartets of a class share one equation, and the
// pseudolikelihood only needs the sum of their observed CFs.
//
// Building the classes costs a few table lookups per quartet, plus one
// equation per class. Equations are found again in the classes of the
// previous network (`reuse`) when the part of the network above the four
// ports is unchanged, whatever the move that was applied.
#pragma once

#include <array>
#include <cstdint>
#include <memory>
#include <unordered_map>
#include <vector>

#include "equations.hpp"
#include "likelihood.hpp"
#include "network.hpp"
#include "optimize.hpp"

namespace csnaq {

// The quartets used by a search: all of them, or a sample.
struct QuartetList {
    int ntaxa = 0;
    bool all = false;                      // every quartet, in rank order
    std::vector<std::array<int, 4>> taxa;  // otherwise: the sorted taxon indices of each quartet
    size_t size() const { return all ? static_cast<size_t>(nchoose4(ntaxa)) : taxa.size(); }
    static QuartetList allquartets(int ntaxa);
    static QuartetList fromranks(int ntaxa, const std::vector<long long>& ranks);
};

// Identifies the part of a network that determines a class's equation: the
// four ports, in order, and every edge above them (128-bit hash).
struct ClassKey {
    uint64_t a = 0, b = 0;
    bool operator==(const ClassKey& o) const { return a == o.a && b == o.b; }
};
struct ClassKeyHash {
    size_t operator()(const ClassKey& k) const { return static_cast<size_t>(k.a ^ (k.b * 0x9e3779b97f4a7c15ULL)); }
};

struct QuartetClasses {
    int npid = 0;
    // per class
    std::vector<ClassKey> key;
    std::vector<std::shared_ptr<const ClassEqn>> eqn;
    std::vector<std::array<double, 3>> Q;  // summed (weighted) observed CFs, in the class's column order
    // pids of the class's four ports (cut edges), in order: the class's columns
    // are the splits 01|23, 02|13 and 03|12 of these ports
    std::vector<std::array<int, 4>> portpid;
    std::unordered_map<ClassKey, uint32_t, ClassKeyHash> bykey;
    // per quartet (only when requested): its class, and the observed-CF column
    // of each class column (2 bits each)
    std::vector<uint32_t> qclass;
    std::vector<uint8_t> qcols;
    // prepared for the fit
    ClassBatch batch;
    size_t nbuilt = 0;  // equations built (the others were reused)
    size_t size() const { return eqn.size(); }

    // Classes whose equation is not built yet (ClassOptions::defer); see finishclasses.
    struct Pending {
        std::vector<uint32_t> classes;
        std::vector<std::array<NodeId, 4>> ports;
        std::vector<NodeId> top;
        std::vector<int> nodepos, edgepos;
        ParamMap pm;
    };
    std::shared_ptr<Pending> pending;
    bool complete() const { return !pending; }
};

struct ClassOptions {
    const double* weights = nullptr;          // weight of each quartet (default: 1)
    const QuartetClasses* reuse = nullptr;    // classes of a previous network to take equations from
    bool perquartet = false;                  // fill qclass / qcols
    double qlogq = 0.0;                       // sum over quartets of weight * q log q (see obsqlogq)
    int nthreads = 1;
    // Do not build the equations of the classes that are not tree-like (the
    // expensive ones): classupperbound can be used right away, and
    // finishclasses builds them if the network is worth fitting.
    bool defer = false;
};

// Sum over the rows of `obs` of weight * sum_k q_k log q_k.
double obsqlogq(const CFMatrix& obs, const double* weights = nullptr);

// The classes of `net` for the quartets `ql`, whose observed CFs are the rows
// of `obs`. Renumbers the nodes and edges of `net`, as SNaQ does before it
// builds quartet equations.
QuartetClasses buildclasses(Network& net, const QuartetList& ql, const CFMatrix& obs, const ClassOptions& opt);

// Builds the equations left out by ClassOptions::defer. `net` must not have
// changed since buildclasses.
void finishclasses(const Network& net, QuartetClasses& cls, int nthreads = 1);

// An upper bound of the SNaQ score of the network over all parameter values:
// each class gets the parameter values that suit it best, whatever the other
// classes need. No fit of the network can score higher.
double classupperbound(const QuartetClasses& cls);

// A pair of edges between which the data suggest a reticulation.
struct HybridCandidate {
    int pid1, pid2;  // Edge::pid of the two edges
    double weight;
};
// Pairs of cut edges that a new hybrid edge could connect to explain what the
// network `net` (with its current parameter values) does not. For every class
// and every split that pairs two of its ports, the observed minus expected
// CFs of the split are added to the pair; the weight of a pair is its total
// excess divided by the square root of its number of quartets. Pairs with a
// positive weight are returned, sorted by decreasing weight.
std::vector<HybridCandidate> hybridcandidates(const Network& net, const QuartetClasses& cls, double rho);

// The value of every parameter of `net`, indexed by pid (`npid` entries).
std::vector<double> pidparams(const Network& net, int npid);

struct FitInfo {
    bool flipped = false;    // a hybrid node's major and minor edges were exchanged
    bool converged = false;  // the optimizer stopped on a tolerance, not on its evaluation limit
    int evals = 0;           // evaluations of the score
};

// fitnumericalparameters! on the classes of `net`.
double fitclasses(Network& net, const QuartetClasses& cls, double rho, const FitOptions& opt, int nthreads = 1,
                  FitInfo* info = nullptr);

// SNaQ score of `net` at its current parameter values.
double scoreclasses(const Network& net, const QuartetClasses& cls, double rho, int nthreads = 1);

}  // namespace csnaq
