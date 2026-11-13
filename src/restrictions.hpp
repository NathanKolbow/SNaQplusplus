// Topological restrictions on candidate networks (network_properties.jl's
// `restrictionset`).
#pragma once

#include <string>

#include "network.hpp"

namespace csnaq {

struct Restrictions {
    int max_level = -1;  // -1: no restriction
    bool galled_tree = false;
    bool galled_network = false;
    bool rooted_tree_child = false;
    bool weakly_tree_child = false;
    bool strongly_tree_child = false;

    // Whether `net` meets all restrictions (restrictionset semantics:
    // galled_tree takes precedence over galled_network; among tree-child
    // options, rooted > weakly > strongly).
    bool operator()(Network& net) const;
    // describerestrictions(...)
    std::string describe() const;
    bool any() const;
};

}  // namespace csnaq
