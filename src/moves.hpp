// Topological moves used by the SNaQ search (network_moves/* and
// optimization/move_proposals.jl).
#pragma once

#include <array>
#include <string>
#include <vector>

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

// samplemoveproposal(N, hmax, rng)
std::pair<Move, MoveParams> samplemoveproposal(Network& net, int hmax, Rng& rng);

// generatemoveproposal(...): returns Move::none if no untried move remains.
std::pair<Move, MoveParams> generatemoveproposal(Network& Nprime, const std::vector<QuartetData>& N_eqns,
                                                 const std::vector<AttemptedMove>& attempted, int hmax,
                                                 double probQR, const CFMatrix& Q, std::vector<double>& CFdeltas,
                                                 Rng& rng, double rho);

// individual moves (exposed for tests)
bool isvalidaddhybrid(const Network& net, EdgeId from_edge, EdgeId to_edge);
NodeId addhybrid(Network& net, EdgeId from_edge, EdgeId to_edge, Rng& rng);
MoveParams sampleaddhybridparameters(Network& net, Rng& rng);

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

}  // namespace csnaq
