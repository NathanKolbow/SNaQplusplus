// Topological moves used by the SNaQ search (network_moves/* and
// optimization/move_proposals.jl).
#pragma once

#include <array>
#include <string>
#include <vector>

#include "classes.hpp"
#include "equations.hpp"
#include "likelihood.hpp"
#include "network.hpp"

namespace csnaq {

enum class Move {
    none,
    addhybrid,
    rNNI1,
    rNNI2,
    rNNI3,
    rNNI4,
    rSPR,
    retic_origin,
    retic_target,
    retic_origin_local,
    retic_target_local,
    flip_hybrid
};

// Julia symbol name (used in logs and for sorting move tables).
const char* move_name(Move m);

// A move's parameters: nodes and/or edges of the network, in order.
struct MoveParam {
    bool is_edge;
    int id;
};
struct MoveParams {
    bool valid = false;  // false: Julia's `nothing`
    std::vector<MoveParam> p;
};

// What `alreadyattempted` compares: the type, number and (for nodes) name of
// each parameter object, as they stand once the iteration that proposed the
// move is over.
struct AttemptedMove {
    Move move;
    std::vector<MoveParam> types;  // id unused
    std::vector<int> numbers;
    std::vector<std::string> names;
};

AttemptedMove snapshot_attempt(const Network& net, Move move, const MoveParams& params);
bool alreadyattempted(const std::vector<AttemptedMove>& attempted, const Network& net, Move move,
                      const MoveParams& params);

bool can_update_inplace(Move m);

// applymove!(N, move, params)
void applymove(Network& net, Move move, const MoveParams& params, Rng& rng);

// Optional guidance for proposals (not in SNaQ.jl).
struct ProposalGuide {
    // Propose only the addition of a hybrid edge.
    bool forceaddhybrid = false;
    // Pairs of edges between which the data suggest a reticulation
    // (hybridcandidates). Half of the hybrid additions are drawn among them,
    // with probability proportional to their weight; the other half
    // uniformly, as SNaQ.jl does.
    const std::vector<HybridCandidate>* candidates = nullptr;
};

// samplemoveproposal(N, hmax, rng)
std::pair<Move, MoveParams> samplemoveproposal(Network& net, int hmax, Rng& rng,
                                               const ProposalGuide* guide = nullptr);

// generatemoveproposal(...): returns Move::none if no untried move remains.
// `classes` are the quartet classes of the current network (of which Nprime is
// a copy), built with their per-quartet data; they are only used when
// probQR > 0 and may be null otherwise.
std::pair<Move, MoveParams> generatemoveproposal(Network& Nprime, const QuartetClasses* classes,
                                                 const std::vector<AttemptedMove>& attempted, int hmax,
                                                 double probQR, const CFMatrix& Q, std::vector<double>& CFdeltas,
                                                 Rng& rng, double rho, const ProposalGuide* guide = nullptr);

// individual moves (exposed for tests)
bool isvalidaddhybrid(const Network& net, EdgeId from_edge, EdgeId to_edge);
NodeId addhybrid(Network& net, EdgeId from_edge, EdgeId to_edge, Rng& rng);
MoveParams sampleaddhybridparameters(Network& net, Rng& rng, const ProposalGuide* guide = nullptr);

bool isvalidrNNI1(const Network& net, NodeId s, NodeId t, NodeId u, NodeId v);
bool isvalidrNNI2(const Network& net, NodeId s, NodeId t, NodeId u, NodeId v);
bool isvalidrNNI3(const Network& net, NodeId s, NodeId t, NodeId u, NodeId v);
bool isvalidrNNI4(const Network& net, NodeId s, NodeId t, NodeId u, NodeId v);
std::vector<std::array<NodeId, 4>> allvalidrNNInodes(const Network& net, int type);
int performrNNI(Network& net, NodeId s, NodeId t, NodeId u, NodeId v, int type, Rng& rng);
MoveParams samplerNNIparameters(const Network& net, int type, Rng& rng);

bool isvalidrSPR(const Network& net, NodeId w, NodeId x, NodeId y, NodeId z, NodeId xprime, NodeId yprime);
void performrSPR(Network& net, NodeId w, NodeId x, NodeId y, NodeId z, NodeId xprime, NodeId yprime);
MoveParams samplerSPRparameters(const Network& net, Rng& rng);

bool isvalidmovereticulateorigin(const Network& net, NodeId hybrid, EdgeId new_origin);
bool isvalidmovereticulatetarget(const Network& net, NodeId hybrid, EdgeId new_target);
bool movereticulateorigin(Network& net, NodeId hybrid, EdgeId new_origin);
bool movereticulatetarget(Network& net, NodeId hybrid, EdgeId new_target);
MoveParams samplemovereticulateoriginparameters(const Network& net, Rng& rng);
MoveParams samplemovereticulatetargetparameters(const Network& net, Rng& rng);
MoveParams samplemovereticulateoriginlocalparameters(const Network& net, Rng& rng);
MoveParams samplemovereticulatetargetlocalparameters(const Network& net, Rng& rng);

MoveParams samplefliphybridparameters(const Network& net, Rng& rng);

// Every local move of `net`, in an order that only depends on the network:
// the NNIs (rNNI1, and rNNI4 if it has reticulations) and, for each
// reticulation, the local moves of its origin and of its target and the change
// of its direction. They are the moves that SnaqOptions::sweep evaluates.
std::vector<std::pair<Move, MoveParams>> localmoves(const Network& net);

}  // namespace csnaq
