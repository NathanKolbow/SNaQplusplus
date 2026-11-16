#include "signature.hpp"

#include <algorithm>
#include <string>
#include <vector>

namespace snaqpp {

namespace {
inline uint64_t mix(uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}
inline uint64_t hashstring(const std::string& s) {
    uint64_t h = 1469598103934665603ULL;
    for (unsigned char c : s) {
        h ^= c;
        h *= 1099511628211ULL;
    }
    return mix(h);
}
}  // namespace

uint64_t topologysignature(const Network& net, bool* discrete) {
    const size_t n = net.node.size();
    std::vector<int> pos(net.npool.size(), -1);
    for (size_t j = 0; j < n; ++j) pos[static_cast<size_t>(net.node[j])] = static_cast<int>(j);
    // neighbors of each node, with the kind of the connecting edge seen from the node:
    // 0 tree edge, 1 hybrid edge into the node, 2 hybrid edge out of the node
    struct Nb {
        int other;
        uint64_t kind;
    };
    std::vector<Nb> nb;
    std::vector<uint32_t> off(n + 1, 0);
    std::vector<uint64_t> color(n), next(n), sorted(n);
    for (size_t j = 0; j < n; ++j) {
        const Node& nd = net.N(net.node[j]);
        for (EdgeId e : nd.edge) {
            const Edge& ed = net.E(e);
            const NodeId o = net.getOtherNode(e, net.node[j]);
            const uint64_t kind = !ed.hybrid ? 0 : (net.getchild(e) == net.node[j] ? 1 : 2);
            nb.push_back({pos[static_cast<size_t>(o)], kind});
        }
        off[j + 1] = static_cast<uint32_t>(nb.size());
        color[j] = nd.leaf ? hashstring(nd.name) : mix(nd.hybrid ? 0x48ULL : 0x54ULL);
    }
    auto ndistinct = [&]() {
        sorted = color;
        std::sort(sorted.begin(), sorted.end());
        return static_cast<size_t>(std::unique(sorted.begin(), sorted.end()) - sorted.begin());
    };
    size_t count = ndistinct();
    for (size_t round = 0; round < n + 1 && count < n; ++round) {
        for (size_t j = 0; j < n; ++j) {
            uint64_t acc = 0;
            for (uint32_t k = off[j]; k < off[j + 1]; ++k)
                acc += mix(color[static_cast<size_t>(nb[k].other)] ^ (0x9e3779b97f4a7c15ULL * (nb[k].kind + 1)));
            next[j] = mix(color[j] ^ mix(acc));
        }
        color.swap(next);
        const size_t c = ndistinct();
        if (c == count) break;  // stable: no further refinement
        count = c;
    }
    // one more round so that every node's color reflects its neighbors' final colors
    for (size_t j = 0; j < n; ++j) {
        uint64_t acc = 0;
        for (uint32_t k = off[j]; k < off[j + 1]; ++k)
            acc += mix(color[static_cast<size_t>(nb[k].other)] ^ (0x9e3779b97f4a7c15ULL * (nb[k].kind + 1)));
        next[j] = mix(color[j] ^ mix(acc));
    }
    std::sort(next.begin(), next.end());
    if (discrete) *discrete = count == n;
    uint64_t h = mix(static_cast<uint64_t>(n) * 1315423911ULL + static_cast<uint64_t>(net.numhybrids));
    for (uint64_t c : next) h = mix(h ^ c);
    return h;
}

}  // namespace snaqpp
