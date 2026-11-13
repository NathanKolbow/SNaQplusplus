#include "restrictions.hpp"

#include <tuple>
#include <vector>

#include "network_manip.hpp"
#include "snaq_network.hpp"

namespace csnaq {

bool Restrictions::any() const {
    return max_level >= 0 || galled_tree || galled_network || rooted_tree_child || weakly_tree_child ||
           strongly_tree_child;
}

bool Restrictions::operator()(Network& net) const {
    if (max_level >= 0 && !(getnetworklevel(net) <= max_level)) return false;
    if (galled_tree) {
        if (!(getnetworklevel(net) <= 1)) return false;
    } else if (galled_network) {
        if (!isgallednetwork(net)) return false;
    }
    if (rooted_tree_child) {
        if (!std::get<0>(istreechild(net))) return false;
    } else if (weakly_tree_child) {
        if (!std::get<1>(istreechild(net))) return false;
    } else if (strongly_tree_child) {
        if (!std::get<2>(istreechild(net))) return false;
    }
    return true;
}

std::string Restrictions::describe() const {
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
