#include "restrictions.hpp"

#include <tuple>
#include <vector>

#include "network_manip.hpp"
#include "snaq_network.hpp"

namespace csnaq {

bool Restrictions::any() const {
    return static_cast<bool>(custom) || max_level >= 0 || galled_tree || galled_network || rooted_tree_child || weakly_tree_child ||
           strongly_tree_child;
}

namespace {

// A copy of `net` whose edges are directed away from its root, with
// Edge::containroot (whether the root can be on the edge) computed again, as
// when the network is written (writenewick) and read. If the current root is
// not admissible, the first admissible one is taken, as writenewick does.
// Returns false if the network has no admissible root.
bool directedcopy(const Network& net, Network& c) {
    c = net.deepcopynetwork();
    const int rootsaved = c.rooti;
    try {
        directedges(c);
        return true;
    } catch (const RootMismatch&) {
    }
    for (EdgeId e : c.edge) {
        c.rooti = c.index_of_node(c.getparent(e));
        try {
            directedges(c);
            return true;
        } catch (const RootMismatch&) {
        }
    }
    c.rooti = rootsaved;
    return false;
}

}  // namespace

bool Restrictions::operator()(Network& net) const {
    if (custom) return custom(net);
    if (!any()) return true;
    // The tree-child conditions read the direction of the edges and
    // Edge::containroot, which the moves of the search do not keep up to date
    // for every edge (SNaQ.jl checks the network as the move left it). The
    // network is checked as it will be read from the output files.
    Network c;
    if (!directedcopy(net, c)) return false;
    if (max_level >= 0 && !(getnetworklevel(c) <= max_level)) return false;
    if (galled_tree) {
        if (!(getnetworklevel(c) <= 1)) return false;
    } else if (galled_network) {
        if (!isgallednetwork(c)) return false;
    }
    if (rooted_tree_child) {
        if (!std::get<0>(istreechild(c))) return false;
    } else if (weakly_tree_child) {
        if (!std::get<1>(istreechild(c))) return false;
    } else if (strongly_tree_child) {
        if (!std::get<2>(istreechild(c))) return false;
    }
    return true;
}

std::string Restrictions::describe() const {
    if (custom) return "user-supplied function `" + custom_name + "`";
    std::vector<std::string> parts;
    if (max_level >= 0) parts.push_back("networks of level ≤ " + std::to_string(max_level));
    if (galled_tree)
        parts.push_back("galled trees (level-1 networks)");
    else if (galled_network)
        parts.push_back("galled networks");
    if (rooted_tree_child)
        parts.push_back("rooted tree-child networks");
    else if (weakly_tree_child)
        parts.push_back("weakly tree-child networks");
    else if (strongly_tree_child)
        parts.push_back("strongly tree-child networks");
    if (parts.empty()) return "none (any network with at most hmax reticulations)";
    std::string out = parts[0];
    for (size_t k = 1; k < parts.size(); ++k) out += " and " + parts[k];
    return out;
}

}  // namespace csnaq
