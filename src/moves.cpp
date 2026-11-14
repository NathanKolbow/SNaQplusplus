#include "moves.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>

#include "network_manip.hpp"
#include "snaq_network.hpp"

namespace csnaq {

const char* move_name(Move m) {
    switch (m) {
        case Move::none: return "none";
        case Move::addhybrid: return "addhybrid";
        case Move::rNNI1: return "rNNI1";
        case Move::rNNI2: return "rNNI2";
        case Move::rNNI3: return "rNNI3";
        case Move::rNNI4: return "rNNI4";
        case Move::rSPR: return "rSPR";
        case Move::retic_origin: return "retic_origin";
        case Move::retic_target: return "retic_target";
        case Move::retic_origin_local: return "retic_origin_local";
        case Move::retic_target_local: return "retic_target_local";
        case Move::flip_hybrid: return "flip_hybrid";
    }
    return "?";
}

bool can_update_inplace(Move m) { return m == Move::rNNI1 || m == Move::rNNI2; }

static MoveParams nodes_params(std::initializer_list<NodeId> ns) {
    MoveParams mp;
    mp.valid = true;
    for (NodeId n : ns) mp.p.push_back({false, n});
    return mp;
}

AttemptedMove snapshot_attempt(const Network& net, Move move, const MoveParams& params) {
    AttemptedMove a;
    a.move = move;
    for (const MoveParam& mp : params.p) {
        a.types.push_back(mp);
        if (mp.is_edge) {
            a.numbers.push_back(net.E(mp.id).number);
            a.names.emplace_back();
        } else {
            a.numbers.push_back(net.N(mp.id).number);
            a.names.push_back(net.N(mp.id).name);
        }
    }
    return a;
}

bool alreadyattempted(const std::vector<AttemptedMove>& attempted, const Network& net, Move move,
                      const MoveParams& params) {
    for (const AttemptedMove& a : attempted) {
        if (a.move != move) continue;
        bool all = true;
        size_t n = std::min(a.types.size(), params.p.size());
        for (size_t k = 0; k < n; ++k) {
            const MoveParam& p = params.p[k];
            if (a.types[k].is_edge != p.is_edge) {
                all = false;
                break;
            }
            int num = p.is_edge ? net.E(p.id).number : net.N(p.id).number;
            if (a.numbers[k] != num || (!p.is_edge && a.names[k] != net.N(p.id).name)) {
                all = false;
                break;
            }
        }
        if (all) return true;
    }
    return false;
}

// ------------------------------------------------------------------------
// add / remove reticulations (add_remove_retic.jl)
// ------------------------------------------------------------------------

bool isvalidaddhybrid(const Network& net, EdgeId from_edge, EdgeId to_edge) {
    if (directionalconflict(net, net.getparent(from_edge), to_edge, true)) return false;
    if (net.getparent(from_edge) == net.getparent(to_edge) && net.isrootof(net.getparent(from_edge))) return false;
    return true;
}

static std::string rand_name_suffix(Rng& rng, long long mod, long long add) {
    return std::to_string(std::llabs(rng.rand_int64() % mod) + add);
}

NodeId addhybrid(Network& net, EdgeId from_edge, EdgeId to_edge, Rng& rng) {
    if (!isvalidaddhybrid(net, from_edge, to_edge)) throw SnaqError("Invalid parameters given to `addhybrid!`");
    auto [newH, newHE] = addhybridedge(net, from_edge, to_edge, true, 0.0, 0.25);
    net.N(net.getparent(from_edge)).name = "i" + rand_name_suffix(rng, 100000, 1000);
    net.N(newH).name = "H" + rand_name_suffix(rng, 100000, 1000);
    // update_yz!(newHE)
    double y = std::exp(-net.E(newHE).length);
    net.E(newHE).y = y;
    net.E(newHE).z = 1.0 - y;
    return newH;
}

MoveParams sampleaddhybridparameters(Network& net, Rng& rng, const ProposalGuide* guide) {
    if (guide && guide->candidates && !guide->candidates->empty() && rng.rand() < 0.5) {
        const std::vector<HybridCandidate>& cands = *guide->candidates;
        std::vector<double> w(cands.size());
        for (size_t k = 0; k < cands.size(); ++k) w[k] = cands[k].weight;
        for (int niter = 0; niter < 20; ++niter) {
            const HybridCandidate& c = cands[rng.weighted_index(w)];
            EdgeId e1 = NONE, e2 = NONE;
            for (EdgeId e : net.edge) {
                if (net.E(e).pid == c.pid1) e1 = e;
                if (net.E(e).pid == c.pid2) e2 = e;
            }
            if (e1 == NONE || e2 == NONE) continue;
            if (rng.rand() < 0.5) std::swap(e1, e2);  // the direction is not known
            MoveParams mp;
            mp.valid = true;
            if (isvalidaddhybrid(net, e1, e2)) {
                mp.p = {{true, e1}, {true, e2}};
                return mp;
            }
            if (isvalidaddhybrid(net, e2, e1)) {
                mp.p = {{true, e2}, {true, e1}};
                return mp;
            }
        }
    }
    for (int niter = 0; niter < 1000; ++niter) {
        auto pick = rng.sample_without_replacement(static_cast<long long>(net.edge.size()), 2);
        EdgeId e1 = net.edge[static_cast<size_t>(pick[0])], e2 = net.edge[static_cast<size_t>(pick[1])];
        MoveParams mp;
        mp.valid = true;
        if (isvalidaddhybrid(net, e1, e2)) {
            mp.p = {{true, e1}, {true, e2}};
            return mp;
        }
        if (isvalidaddhybrid(net, e2, e1)) {
            mp.p = {{true, e2}, {true, e1}};
            return mp;
        }
    }
    return {};
}

// ------------------------------------------------------------------------
// rNNI moves (rNNI_moves.jl, rNNI_validity.jl)
// ------------------------------------------------------------------------

template <class V>
static bool contains_node(const V& v, NodeId n) {
    return std::find(v.begin(), v.end(), n) != v.end();
}

static EdgeId edge_with(const Network& net, NodeId a, NodeId b) {
    for (EdgeId e : net.N(a).edge)
        if (net.E(e).has_node(b)) return e;
    throw SnaqError("expected nodes to be adjacent");
}

static void replace_edge_node(Network& net, EdgeId e, NodeId from, NodeId to) {
    Edge& ed = net.E(e);
    for (int k = 0; k < ed.nnode; ++k)
        if (ed.nd[k] == from) ed.nd[k] = to;
}

static bool distinct4(NodeId s, NodeId t, NodeId u, NodeId v) {
    return !(s == t || s == u || s == v || t == u || t == v || u == v);
}

bool isvalidrNNI1(const Network& net, NodeId s, NodeId t, NodeId u, NodeId v) {
    if (!distinct4(s, t, u, v)) return false;
    if (net.N(s).hybrid || net.N(t).hybrid || net.N(u).hybrid || net.N(v).hybrid) {
        auto cu = net.getchildren(u);
        if (!(contains_node(cu, s) && contains_node(cu, v))) return false;
        if (isdescendantof(net, s, v)) return false;
        auto cv = net.getchildren(v);
        if (!(!contains_node(cv, s) && !contains_node(net.getchildren(s), v))) return false;
    }
    return areconnected(net, u, v) && areconnected(net, s, u) && areconnected(net, t, v);
}

bool isvalidrNNI2(const Network& net, NodeId s, NodeId t, NodeId u, NodeId v) {
    if (!distinct4(s, t, u, v)) return false;
    if (!contains_node(net.getchildren(u), v)) return false;
    if (!contains_node(net.getchildren(t), v)) return false;
    if (!contains_node(net.getchildren(s), u)) return false;
    if (!net.N(v).hybrid) return false;
    if (isdescendantof(net, t, u)) return false;
    return true;
}

bool isvalidrNNI3(const Network& net, NodeId s, NodeId t, NodeId u, NodeId v) {
    if (!distinct4(s, t, u, v)) return false;
    if (net.N(v).hybrid) return false;
    if (!contains_node(net.getchildren(v), t)) return false;
    if (!contains_node(net.getchildren(s), u)) return false;
    if (!contains_node(net.getchildren(u), v)) return false;
    if (!net.N(u).hybrid) return false;
    auto pu = net.getparents(u);
    if (pu.size() < 2) throw SnaqError("BoundsError: hybrid node with fewer than 2 parents");
    if (pu[0] == pu[1]) return false;
    return true;
}

bool isvalidrNNI4(const Network& net, NodeId s, NodeId t, NodeId u, NodeId v) {
    if (!distinct4(s, t, u, v)) return false;
    if (!contains_node(net.getchildren(u), v)) return false;
    if (!contains_node(net.getchildren(t), v)) return false;
    if (!contains_node(net.getchildren(u), s)) return false;
    if (net.getchildren(u).size() < 2) return false;
    if (!net.N(v).hybrid) return false;
    if (isdescendantof(net, u, t)) return false;
    if (isdescendantof(net, t, s)) return false;
    return true;
}

static std::vector<std::array<NodeId, 4>> allvalidrNNI1nodes(const Network& net) {
    std::vector<NodeId> nonleaves;
    for (NodeId n : net.node)
        if (!net.N(n).leaf && !net.N(n).hybrid) nonleaves.push_back(n);
    std::vector<std::array<NodeId, 4>> out;
    for (size_t i = 0; i + 1 < nonleaves.size(); ++i) {
        NodeId u = nonleaves[i];
        for (size_t j = i + 1; j < nonleaves.size(); ++j) {
            NodeId v = nonleaves[j];
            if (!areconnected(net, u, v)) continue;
            for (EdgeId se : net.N(u).edge) {
                NodeId s = net.getOtherNode(se, u);
                if (s == v || net.N(s).hybrid) continue;
                for (EdgeId te : net.N(v).edge) {
                    NodeId t = net.getOtherNode(te, v);
                    if (t == u || net.N(t).hybrid) continue;
                    out.push_back({s, t, u, v});
                }
            }
        }
    }
    return out;
}

static std::vector<std::array<NodeId, 4>> allvalidrNNI2nodes(const Network& net) {
    std::vector<std::array<NodeId, 4>> out;
    for (NodeId v : net.node) {
        if (!net.N(v).hybrid) continue;
        auto parents = net.getparents(v);
        if (parents[0] == parents[1]) continue;
        auto p1p = net.getparents(parents[0]);
        auto p2p = net.getparents(parents[1]);
        if (!p1p.empty() && !isdescendantof(net, parents[1], parents[0]))
            for (NodeId gpa : p1p) {
                if (isdescendantof(net, gpa, v)) continue;
                if (gpa == parents[1] || parents[0] == v) continue;
                out.push_back({gpa, parents[1], parents[0], v});
            }
        if (!p2p.empty() && !isdescendantof(net, parents[0], parents[1]))
            for (NodeId gpa : p2p) {
                if (isdescendantof(net, gpa, v)) continue;
                if (gpa == parents[0] || parents[1] == v) continue;
                out.push_back({gpa, parents[0], parents[1], v});
            }
    }
    for (auto& c : out)
        if (!isvalidrNNI2(net, c[0], c[1], c[2], c[3])) throw SnaqError("Some generated combos weren't valid...");
    return out;
}

static std::vector<std::array<NodeId, 4>> allvalidrNNI3nodes(const Network& net) {
    std::vector<std::array<NodeId, 4>> out;
    for (NodeId u : net.node) {
        if (!net.N(u).hybrid || net.N(net.getchild_node(u)).leaf) continue;
        NodeId v = net.getchild_node(u);
        if (net.N(v).hybrid) continue;
        auto ss = net.getparents(u);
        if (ss[0] == ss[1]) continue;
        for (NodeId s : ss)
            for (NodeId t : net.getchildren(v)) out.push_back({s, t, u, v});
    }
    for (auto& c : out)
        if (!isvalidrNNI3(net, c[0], c[1], c[2], c[3])) throw SnaqError("Some generated combos weren't valid...");
    return out;
}

static std::vector<std::array<NodeId, 4>> allvalidrNNI4nodes(const Network& net) {
    std::vector<std::array<NodeId, 4>> out;
    for (NodeId v : net.node) {
        if (!net.N(v).hybrid) continue;
        auto parents = net.getparents(v);
        auto p1c = net.getchildren(parents[0]);
        auto p2c = net.getchildren(parents[1]);
        if (p1c.size() == 2 && !isdescendantof(net, parents[0], parents[1]))
            for (NodeId s : p1c) {
                if (s == parents[1] || s == v) continue;
                if (isdescendantof(net, parents[1], s)) continue;
                out.push_back({s, parents[1], parents[0], v});
            }
        if (p2c.size() == 2 && !isdescendantof(net, parents[1], parents[0]))
            for (NodeId s : p2c) {
                if (s == parents[0] || s == v) continue;
                if (isdescendantof(net, parents[0], s)) continue;
                out.push_back({s, parents[0], parents[1], v});
            }
    }
    for (auto& c : out)
        if (!isvalidrNNI4(net, c[0], c[1], c[2], c[3])) throw SnaqError("Some generated combos weren't valid...");
    return out;
}

std::vector<std::array<NodeId, 4>> allvalidrNNInodes(const Network& net, int type) {
    if (type == 1) return allvalidrNNI1nodes(net);
    if (type == 2) return allvalidrNNI2nodes(net);
    if (type == 3) return allvalidrNNI3nodes(net);
    return allvalidrNNI4nodes(net);
}

static std::string hybridize_name(const std::string& name, Rng& rng) {
    if (name.empty()) return "H" + std::to_string(std::llabs(rng.rand_int64() % 10000) + 100);
    if (name[0] == 'H') return name;
    if (name[0] == 'i') return "H" + name.substr(1);
    return "H" + name;
}

static std::string treeify_name(const std::string& name) {
    // QUIRK: SNaQ.jl indexes name[1] here, which errors if the name is empty.
    if (name.empty()) throw SnaqError("BoundsError: attempt to access empty String at index [1]");
    return name[0] == 'H' ? "i" + name.substr(1) : name;
}

static void performrNNI1(Network& net, NodeId s, NodeId t, NodeId u, NodeId v) {
    if (!isvalidrNNI1(net, s, t, u, v)) throw SnaqError("Topological conditions for rNNI(1) not met.");
    EdgeId edge_su = edge_with(net, s, u);
    EdgeId edge_tv = edge_with(net, t, v);
    replace_edge_node(net, edge_su, s, t);
    replace_edge_node(net, edge_tv, t, s);
    net.N(s).edge.replace(edge_su, edge_tv);
    net.N(t).edge.replace(edge_tv, edge_su);
    if (net.N(s).hybrid && net.getchild(edge_tv) != s) std::swap(net.E(edge_tv).nd[0], net.E(edge_tv).nd[1]);
    if (net.N(t).hybrid && net.getchild(edge_su) != t) std::swap(net.E(edge_su).nd[0], net.E(edge_su).nd[1]);
    swapedgehybridinfo(net, edge_su, edge_tv);
    directedges(net);
}

static void performrNNI2(Network& net, NodeId s, NodeId t, NodeId u, NodeId v) {
    if (!isvalidrNNI2(net, s, t, u, v)) throw SnaqError("Topological conditions for rNNI(2) not met.");
    EdgeId edge_su = edge_with(net, s, u);
    EdgeId edge_tv = edge_with(net, t, v);
    EdgeId edge_uv = edge_with(net, u, v);
    if (!(net.E(edge_uv).hybrid && net.E(edge_tv).hybrid)) throw SnaqError("Topological conditions for move not met");
    replace_edge_node(net, edge_su, u, v);
    replace_edge_node(net, edge_tv, v, u);
    net.N(u).edge.replace(edge_su, edge_tv);
    net.N(v).edge.replace(edge_tv, edge_su);
    swapedgehybridinfo(net, edge_su, edge_tv);
}

static EdgeId edge_avoiding(const Network& net, NodeId u, NodeId a, NodeId b) {
    for (EdgeId e : net.N(u).edge)
        if (!net.E(e).has_node(a) && !net.E(e).has_node(b)) return e;
    throw SnaqError("no edge found");
}

static void performrNNI3(Network& net, NodeId s, NodeId t, NodeId u, NodeId v, Rng& rng) {
    if (!isvalidrNNI3(net, s, t, u, v)) throw SnaqError("Topological conditions for rNNI(3) not met.");
    EdgeId edge_su = edge_with(net, s, u);
    EdgeId edge_tv = edge_with(net, t, v);
    EdgeId edge_uv = edge_with(net, u, v);
    EdgeId edge_u = edge_avoiding(net, u, s, v);
    if (!net.E(edge_su).hybrid) throw SnaqError("Topological conditions for move not met");
    replace_edge_node(net, edge_su, u, v);
    replace_edge_node(net, edge_tv, v, u);
    net.N(u).edge.replace(edge_su, edge_tv);
    net.N(v).edge.replace(edge_tv, edge_su);
    swapedgehybridinfo(net, edge_u, edge_uv);
    net.N(u).hybrid = false;
    net.N(v).hybrid = true;
    net.N(v).name = hybridize_name(net.N(v).name, rng);
    net.N(u).name = treeify_name(net.N(u).name);
    std::replace(net.hybrid.begin(), net.hybrid.end(), u, v);
}

static void performrNNI4(Network& net, NodeId s, NodeId t, NodeId u, NodeId v, Rng& rng) {
    if (!isvalidrNNI4(net, s, t, u, v)) throw SnaqError("Topological conditions for rNNI(4) not met.");
    EdgeId edge_su = edge_with(net, s, u);
    EdgeId edge_tv = edge_with(net, t, v);
    EdgeId edge_uv = edge_with(net, u, v);
    EdgeId edge_u = edge_avoiding(net, u, s, v);
    replace_edge_node(net, edge_su, u, v);
    replace_edge_node(net, edge_tv, v, u);
    net.N(u).edge.replace(edge_su, edge_tv);
    net.N(v).edge.replace(edge_tv, edge_su);
    swapedgehybridinfo(net, edge_u, edge_uv);
    net.N(u).hybrid = true;
    net.N(v).hybrid = false;
    net.N(u).name = hybridize_name(net.N(u).name, rng);
    net.N(v).name = treeify_name(net.N(v).name);
    std::replace(net.hybrid.begin(), net.hybrid.end(), v, u);
}

int performrNNI(Network& net, NodeId s, NodeId t, NodeId u, NodeId v, int type, Rng& rng) {
    switch (type) {
        case 1: performrNNI1(net, s, t, u, v); return 1;
        case 2: performrNNI2(net, s, t, u, v); return 2;
        case 3: performrNNI3(net, s, t, u, v, rng); return 3;
        case 4: performrNNI4(net, s, t, u, v, rng); return 4;
    }
    return 0;
}

MoveParams samplerNNIparameters(const Network& net, int type, Rng& rng) {
    if (type < 1 || type > 4) throw SnaqError("`type` must be 1, 2, 3, or 4.");
    auto valid = allvalidrNNInodes(net, type);
    if (valid.empty()) return {};
    const auto& c = valid[rng.index(valid.size())];
    return nodes_params({c[0], c[1], c[2], c[3]});
}

// ------------------------------------------------------------------------
// rSPR (rSPR_moves.jl, rSPR_validity.jl)
// ------------------------------------------------------------------------

bool isvalidrSPR(const Network& net, NodeId w, NodeId x, NodeId y, NodeId z, NodeId xprime, NodeId yprime) {
    if (net.N(z).hybrid) {
        if (!contains_node(net.getchildren(w), z)) return false;
    } else {
        if (!contains_node(net.getchildren(z), w)) return false;
    }
    if (!contains_node(net.getchildren(xprime), yprime)) return false;
    std::vector<NodeId> all{w, x, y, z, xprime, yprime};
    std::sort(all.begin(), all.end());
    size_t nuniq = static_cast<size_t>(std::unique(all.begin(), all.end()) - all.begin());
    if (nuniq != 6) {
        if (!(nuniq == 5 && (y == xprime || x == xprime))) return false;
    }
    if (net.N(z).hybrid) {
        if (isdescendantof(net, w, yprime)) return false;
    } else {
        if (isdescendantof(net, xprime, w)) return false;
    }
    return true;
}

void performrSPR(Network& net, NodeId w, NodeId x, NodeId y, NodeId z, NodeId xprime, NodeId yprime) {
    if (!isvalidrSPR(net, w, x, y, z, xprime, yprime)) throw SnaqError("Topological conditions for rSPR not met.");
    const int which_case = net.N(z).hybrid ? 2 : 1;
    EdgeId edge_xprime_yprime = edge_with(net, xprime, yprime);
    auto [zprime, edge_xprime_zprime] = breakedge(net, edge_xprime_yprime);
    net.N(zprime).name = net.N(z).name;
    net.N(zprime).number = net.N(z).number;
    EdgeId edge_zw = edge_with(net, z, w);
    replace_edge_node(net, edge_zw, z, zprime);
    net.N(z).edge.erase_at(static_cast<size_t>(net.N(z).edge.find(edge_zw)));
    net.N(zprime).edge.push_back(edge_zw);
    if (which_case == 2) {
        EdgeId edge_xz = edge_with(net, z, x);
        swapedgehybridinfo(net, edge_xprime_zprime, edge_xz);
        net.N(zprime).hybrid = true;
        net.pushHybrid(zprime);
    }
    if (net.N(z).hybrid) {
        net.removeHybrid(z);
        net.N(z).hybrid = false;
        for (EdgeId e : net.N(z).edge) {
            net.E(e).hybrid = false;
            net.E(e).gamma = -1.0;
        }
    }
    int z_idx = net.index_of_node(z);
    EdgeId fuseedge = fuseedgesat(net, z_idx);
    NodeId fc = net.getchild(fuseedge);
    if (net.N(fc).hybrid) {
        EdgeId other = NONE;
        for (EdgeId e : net.N(fc).edge)
            if (net.E(e).hybrid && e != fuseedge && net.getchild(e) == fc) {
                other = e;
                break;
            }
        if (other == NONE) throw SnaqError("BoundsError: no partner hybrid edge after rSPR fuse");
        net.E(fuseedge).hybrid = true;
        net.E(fuseedge).gamma = 1 - net.E(other).gamma;
        net.E(fuseedge).ismajor = net.E(fuseedge).gamma > 0.5;
        net.E(other).ismajor = net.E(other).gamma >= 0.5;
    }
    if (std::find(net.hybrid.begin(), net.hybrid.end(), z) != net.hybrid.end())
        throw SnaqError("z still in net.hybrid after being fused??");
}

MoveParams samplerSPRparameters(const Network& net, Rng& rng) {
    long attempts = 0;
    while (true) {
        ++attempts;
        if (attempts >= 100000) return {};
        NodeId xprime = net.node[rng.index(net.node.size())];
        if (net.N(xprime).leaf) continue;
        auto xc = net.getchildren(xprime);
        NodeId yprime = (xc.size() == 1) ? net.getchild_node(xprime) : xc[rng.index(xc.size())];
        NodeId z = net.node[rng.index(net.node.size())];
        if (net.N(z).leaf) continue;
        if (net.N(z).hybrid) {
            NodeId w = net.getparent(net.getparentedgeminor(z));
            NodeId x = net.getparent(net.getparentedge(z));
            NodeId y = net.getchild_node(z);
            if (isvalidrSPR(net, w, x, y, z, xprime, yprime)) return nodes_params({w, x, y, z, xprime, yprime});
        } else {
            auto zc = net.getchildren(z);
            if (zc.size() != 2) continue;
            NodeId w = zc[0], y = zc[1];
            NodeId x = net.getparent_node(z);
            if (isvalidrSPR(net, w, x, y, z, xprime, yprime)) return nodes_params({w, x, y, z, xprime, yprime});
        }
    }
}

// ------------------------------------------------------------------------
// moving reticulation origins / targets (move_origin_target.jl)
// ------------------------------------------------------------------------

// The (x, y) pair for moving the origin of `hybrid`; returns false if, at the
// root, there are not exactly 2 other neighbors.
static bool origin_xy(const Network& net, NodeId hybrid, NodeId z, NodeId& x, NodeId& y) {
    if (!net.getparents(z).empty()) {
        auto zc = net.getchildren(z);
        y = (zc.size() == 1 || zc[0] != hybrid) ? zc[0] : zc[1];
        x = net.getparent_node(z);
        return true;
    }
    std::vector<NodeId> others;
    for (EdgeId e : net.N(z).edge) {
        NodeId o = net.getOtherNode(e, z);
        if (o != hybrid) others.push_back(o);
    }
    if (others.size() != 2) return false;
    x = others[0];
    y = others[1];
    return true;
}

bool isvalidmovereticulateorigin(const Network& net, NodeId hybrid, EdgeId new_origin) {
    if (!net.N(hybrid).hybrid) return false;
    NodeId z = net.getparent(net.getparentedgeminor(hybrid));
    if (net.N(z).edge.size() != 3)
        throw SnaqError("Node with only " + std::to_string(net.N(z).edge.size()) + " attached edges...");
    NodeId xprime = net.getparent(new_origin), yprime = net.getchild(new_origin);
    NodeId x, y;
    if (!origin_xy(net, hybrid, z, x, y)) return false;
    return isvalidrSPR(net, hybrid, x, y, z, xprime, yprime);
}

bool movereticulateorigin(Network& net, NodeId hybrid, EdgeId new_origin) {
    if (!net.N(hybrid).hybrid) throw SnaqError("`hybrid` is not a hybrid node...");
    NodeId w = hybrid;
    NodeId z = net.getparent(net.getparentedgeminor(hybrid));
    NodeId x, y;
    if (!origin_xy(net, hybrid, z, x, y)) return false;
    NodeId xprime = net.getparent(new_origin), yprime = net.getchild(new_origin);
    if (!isvalidrSPR(net, w, x, y, z, xprime, yprime)) throw SnaqError("Invalid `movereticulateorigin!` move.");
    performrSPR(net, w, x, y, z, xprime, yprime);
    return true;
}

bool isvalidmovereticulatetarget(const Network& net, NodeId hybrid, EdgeId new_target) {
    if (!net.N(hybrid).hybrid) throw SnaqError("`hybrid` is not a hybrid node...");
    NodeId x = net.getparent(net.getparentedge(hybrid));
    NodeId w = net.getparent(net.getparentedgeminor(hybrid));
    NodeId y = net.getchild_node(hybrid);
    return isvalidrSPR(net, w, x, y, hybrid, net.getparent(new_target), net.getchild(new_target));
}

bool movereticulatetarget(Network& net, NodeId hybrid, EdgeId new_target) {
    if (!net.N(hybrid).hybrid) throw SnaqError("`hybrid` is not a hybrid node...");
    NodeId x = net.getparent(net.getparentedge(hybrid));
    NodeId w = net.getparent(net.getparentedgeminor(hybrid));
    NodeId y = net.getchild_node(hybrid);
    NodeId xprime = net.getparent(new_target), yprime = net.getchild(new_target);
    if (!isvalidrSPR(net, w, x, y, hybrid, xprime, yprime)) throw SnaqError("Invalid `movereticulatetarget!` move.");
    performrSPR(net, w, x, y, hybrid, xprime, yprime);
    return true;
}

static MoveParams hyb_edge_params(NodeId h, EdgeId e) {
    MoveParams mp;
    mp.valid = true;
    mp.p = {{false, h}, {true, e}};
    return mp;
}

MoveParams samplemovereticulateoriginparameters(const Network& net, Rng& rng) {
    for (int hi : rng.randperm(static_cast<int>(net.hybrid.size()))) {
        NodeId H = net.hybrid[static_cast<size_t>(hi)];
        for (int ei : rng.randperm(static_cast<int>(net.edge.size()))) {
            EdgeId e = net.edge[static_cast<size_t>(ei)];
            if (isvalidmovereticulateorigin(net, H, e)) return hyb_edge_params(H, e);
        }
    }
    return {};
}

MoveParams samplemovereticulatetargetparameters(const Network& net, Rng& rng) {
    for (int hi : rng.randperm(static_cast<int>(net.hybrid.size()))) {
        NodeId H = net.hybrid[static_cast<size_t>(hi)];
        for (int ei : rng.randperm(static_cast<int>(net.edge.size()))) {
            EdgeId e = net.edge[static_cast<size_t>(ei)];
            if (isvalidmovereticulatetarget(net, H, e)) return hyb_edge_params(H, e);
        }
    }
    return {};
}

MoveParams samplemovereticulateoriginlocalparameters(const Network& net, Rng& rng) {
    if (net.numhybrids <= 0) throw SnaqError("`N` must have at least 1 hybrid node.");
    for (int hi : rng.randperm(net.numhybrids)) {
        NodeId hybrid = net.hybrid[static_cast<size_t>(hi)];
        EdgeId minor_edge = net.getparentedgeminor(hybrid);
        NodeId origin = net.getparent(minor_edge);
        std::vector<EdgeId> cand;
        for (EdgeId e : net.N(origin).edge) {
            if (e == minor_edge) continue;
            NodeId o = net.getOtherNode(e, origin);
            for (EdgeId ne : net.N(o).edge)
                if (ne != e) cand.push_back(ne);
        }
        for (int j : rng.randperm(static_cast<int>(cand.size())))
            if (isvalidmovereticulateorigin(net, hybrid, cand[static_cast<size_t>(j)]))
                return hyb_edge_params(hybrid, cand[static_cast<size_t>(j)]);
    }
    return {};
}

MoveParams samplemovereticulatetargetlocalparameters(const Network& net, Rng& rng) {
    if (net.numhybrids <= 0) throw SnaqError("`N` must have at least 1 hybrid node.");
    for (int hi : rng.randperm(net.numhybrids)) {
        NodeId hybrid = net.hybrid[static_cast<size_t>(hi)];
        NodeId major_parent = net.getparent(net.getparentedge(hybrid));
        std::vector<EdgeId> cand;
        for (EdgeId e : net.N(major_parent).edge)
            if (e != net.getparentedge(hybrid)) cand.push_back(e);
        auto hc = net.getchildren(hybrid);
        if (hc.size() == 1) hc = net.getchildren(hc[0]);
        for (NodeId c : hc) cand.push_back(net.getparentedge(c));
        for (int j : rng.randperm(static_cast<int>(cand.size())))
            if (isvalidmovereticulatetarget(net, hybrid, cand[static_cast<size_t>(j)]))
                return hyb_edge_params(hybrid, cand[static_cast<size_t>(j)]);
    }
    return {};
}

MoveParams samplefliphybridparameters(const Network& net, Rng& rng) {
    for (int j : rng.randperm(net.numhybrids)) {
        Network n0 = net.deepcopynetwork();
        try {
            if (!fliphybrid(n0, n0.hybrid[static_cast<size_t>(j)])) continue;
            return nodes_params({net.hybrid[static_cast<size_t>(j)]});
        } catch (const std::exception&) {
            continue;
        }
    }
    return {};
}

std::vector<std::pair<Move, MoveParams>> localmoves(const Network& net) {
    std::vector<std::pair<Move, MoveParams>> all;
    for (const auto& c : allvalidrNNInodes(net, 1)) all.emplace_back(Move::rNNI1, nodes_params({c[0], c[1], c[2], c[3]}));
    if (net.numhybrids == 0) return all;
    for (const auto& c : allvalidrNNInodes(net, 4)) all.emplace_back(Move::rNNI4, nodes_params({c[0], c[1], c[2], c[3]}));
    for (NodeId hybrid : net.hybrid) {
        // as samplemovereticulateoriginlocalparameters
        const EdgeId minor_edge = net.getparentedgeminor(hybrid);
        const NodeId origin = net.getparent(minor_edge);
        for (EdgeId e : net.N(origin).edge) {
            if (e == minor_edge) continue;
            const NodeId o = net.getOtherNode(e, origin);
            for (EdgeId ne : net.N(o).edge)
                if (ne != e && isvalidmovereticulateorigin(net, hybrid, ne))
                    all.emplace_back(Move::retic_origin_local, hyb_edge_params(hybrid, ne));
        }
        // as samplemovereticulatetargetlocalparameters
        const NodeId major_parent = net.getparent(net.getparentedge(hybrid));
        std::vector<EdgeId> cand;
        for (EdgeId e : net.N(major_parent).edge)
            if (e != net.getparentedge(hybrid)) cand.push_back(e);
        auto hc = net.getchildren(hybrid);
        if (hc.size() == 1) hc = net.getchildren(hc[0]);
        for (NodeId c : hc) cand.push_back(net.getparentedge(c));
        for (EdgeId e : cand)
            if (isvalidmovereticulatetarget(net, hybrid, e))
                all.emplace_back(Move::retic_target_local, hyb_edge_params(hybrid, e));
        // as samplefliphybridparameters
        Network n0 = net.deepcopynetwork();
        bool flips = false;
        try {
            flips = fliphybrid(n0, hybrid);
        } catch (const std::exception&) {
        }
        if (flips) all.emplace_back(Move::flip_hybrid, nodes_params({hybrid}));
    }
    return all;
}

// ------------------------------------------------------------------------
// proposals (move_proposals.jl)
// ------------------------------------------------------------------------

void applymove(Network& net, Move move, const MoveParams& mp, Rng& rng) {
    auto id = [&](size_t k) { return mp.p.at(k).id; };
    switch (move) {
        case Move::addhybrid: addhybrid(net, id(0), id(1), rng); return;
        case Move::rNNI1: performrNNI(net, id(0), id(1), id(2), id(3), 1, rng); return;
        case Move::rNNI2: performrNNI(net, id(0), id(1), id(2), id(3), 2, rng); return;
        case Move::rNNI3: performrNNI(net, id(0), id(1), id(2), id(3), 3, rng); return;
        case Move::rNNI4: performrNNI(net, id(0), id(1), id(2), id(3), 4, rng); return;
        case Move::retic_origin:
        case Move::retic_origin_local: movereticulateorigin(net, id(0), id(1)); return;
        case Move::retic_target:
        case Move::retic_target_local: movereticulatetarget(net, id(0), id(1)); return;
        case Move::rSPR: performrSPR(net, id(0), id(1), id(2), id(3), id(4), id(5)); return;
        case Move::flip_hybrid: fliphybrid(net, id(0)); return;
        case Move::none: break;
    }
    throw SnaqError(std::string("Move \"") + move_name(move) + "\" not recognized.");
}

std::pair<Move, MoveParams> samplemoveproposal(Network& net, int hmax, Rng& rng, const ProposalGuide* guide) {
    if (net.numhybrids < hmax && ((guide && guide->forceaddhybrid) || rng.rand() < 0.05))
        return {Move::addhybrid, sampleaddhybridparameters(net, rng, guide)};
    if (net.numhybrids == 0) {
        double r = rng.rand();
        if (r < 0.7) return {Move::rNNI1, samplerNNIparameters(net, 1, rng)};
        return {Move::rSPR, samplerSPRparameters(net, rng)};
    }
    // rNNI1 15%, rNNI2 0%, rNNI3 0%, rNNI4 10%, rSPR 5%, origin 10%, target 10%,
    // local origin 20%, local target 15%, fliphybrid 15%
    const double probs[10] = {0.15, 0.0, 0.0, 0.1, 0.05, 0.1, 0.1, 0.2, 0.15, 0.15};
    double cum[10];
    double acc = 0.0;
    for (int k = 0; k < 10; ++k) {
        acc += probs[k];
        cum[k] = acc;
    }
    double r = rng.rand();
    if (r <= cum[0]) return {Move::rNNI1, samplerNNIparameters(net, 1, rng)};
    if (r <= cum[1]) return {Move::rNNI2, samplerNNIparameters(net, 2, rng)};
    if (r <= cum[2]) return {Move::rNNI3, samplerNNIparameters(net, 3, rng)};
    if (r <= cum[3]) return {Move::rNNI4, samplerNNIparameters(net, 4, rng)};
    if (r <= cum[4]) return {Move::rSPR, samplerSPRparameters(net, rng)};
    if (r <= cum[5]) return {Move::retic_origin, samplemovereticulateoriginparameters(net, rng)};
    if (r <= cum[6]) return {Move::retic_target, samplemovereticulatetargetparameters(net, rng)};
    if (r <= cum[7]) return {Move::retic_origin_local, samplemovereticulateoriginlocalparameters(net, rng)};
    if (r <= cum[8]) return {Move::retic_target_local, samplemovereticulatetargetlocalparameters(net, rng)};
    return {Move::flip_hybrid, samplefliphybridparameters(net, rng)};
}

static EdgeId sampleprobQRedge(Network& N, const QuartetClasses& cls, const CFMatrix& Q,
                               std::vector<double>& CFdeltas, Rng& rng, double rho) {
    OptInfo info = gatheroptimizationinfo(N, false);
    if (CFdeltas.empty()) {
        // QUIRK: gatherparams(N) renumbers N (gatheroptimizationinfo(N, true)).
        gatheroptimizationinfo(N, true);
        const std::vector<double> params = pidparams(N, cls.npid);
        const double alpha = rhotoalpha(rho);
        std::vector<std::array<double, 3>> ecf(cls.size());
        for (size_t c = 0; c < cls.size(); ++c) ecf[c] = classexpectedCF(*cls.eqn[c], params.data(), cls.npid, alpha);
        for (size_t j = 0; j < cls.qclass.size(); ++j) {
            const auto& e = ecf[cls.qclass[j]];
            const uint8_t cols = cls.qcols[j];
            double d = 0.0;
            for (int k = 0; k < 3; ++k) d += std::abs(e[static_cast<size_t>(k)] - Q(j, (cols >> (2 * k)) & 3));
            CFdeltas.push_back(d);
        }
    }
    size_t iq = rng.weighted_index(CFdeltas);
    std::vector<const RecursiveCFEquation*> queue{&cls.eqn[cls.qclass[iq]]->eqn};
    std::vector<int> edges;
    size_t head = 0;
    while (head < queue.size()) {
        const RecursiveCFEquation* curr = queue[head++];
        for (int p : curr->coal_edges) edges.push_back(p);
        for (const auto& d : curr->divisions) queue.push_back(&d);
    }
    if (!edges.empty()) {
        std::vector<int> uniq;
        for (int p : edges)
            if (std::find(uniq.begin(), uniq.end(), p) == uniq.end()) uniq.push_back(p);
        const int pid = uniq[rng.index(uniq.size())];
        for (EdgeId e : N.edge)
            if (N.E(e).pid == pid) return e;
        throw SnaqError("sampleprobQRedge: the classes do not belong to this network");
    }
    std::vector<int> validkeys;
    for (size_t k = 0; k < info.idx_obj.size(); ++k)
        if (!info.idx_obj[k].is_hybrid_node) validkeys.push_back(static_cast<int>(k));
    return info.idx_obj.at(static_cast<size_t>(validkeys[rng.index(validkeys.size())])).id;
}

std::pair<Move, MoveParams> generatemoveproposal(Network& Nprime, const QuartetClasses* classes,
                                                 const std::vector<AttemptedMove>& attempted, int hmax,
                                                 double probQR, const CFMatrix& Q, std::vector<double>& CFdeltas,
                                                 Rng& rng, double rho, const ProposalGuide* guide) {
    EdgeId required_edge = NONE;
    if (probQR > 0.0 && (!classes || classes->qclass.empty()))
        throw SnaqError("generatemoveproposal: probQR > 0 needs the quartet classes of the current network");
    if (probQR > 0.0 && rng.rand() <= probQR) required_edge = sampleprobQRedge(Nprime, *classes, Q, CFdeltas, rng, rho);
    auto validmove = [&](Move mv, const MoveParams& pars) {
        if (!pars.valid) return false;
        if (alreadyattempted(attempted, Nprime, mv, pars)) return false;
        if (required_edge == NONE) return true;
        for (const MoveParam& p : pars.p)
            if (p.is_edge && p.id == required_edge) return true;
        return false;
    };
    long retries = 0;
    auto mp = samplemoveproposal(Nprime, hmax, rng, guide);
    int attempted_reqedges = 1;
    while (!validmove(mp.first, mp.second)) {
        mp = samplemoveproposal(Nprime, hmax, rng, guide);
        ++retries;
        if (retries >= 1000) {
            if (required_edge == NONE || attempted_reqedges >= 100) return {Move::none, {}};
            required_edge = sampleprobQRedge(Nprime, *classes, Q, CFdeltas, rng, rho);
            ++attempted_reqedges;
            retries = 0;
        }
    }
    return mp;
}

}  // namespace csnaq
