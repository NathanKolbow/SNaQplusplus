#include "equations.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <unordered_map>
#include <unordered_set>

#include "network_manip.hpp"
#include "parallel.hpp"
#include "quartets.hpp"

namespace csnaq {

// ------------------------------------------------------------------------
// ParamMap
// ------------------------------------------------------------------------

void ParamMap::set(int number, int idx) {
    if (v_.empty()) {
        base_ = number;
        v_.assign(1, -1);
    }
    long long k = static_cast<long long>(number) - base_;
    if (k < 0) {
        std::vector<int> nv(static_cast<size_t>(v_.size() - k), -1);
        std::copy(v_.begin(), v_.end(), nv.begin() + (-k));
        v_.swap(nv);
        base_ = number;
        k = 0;
    } else if (k >= static_cast<long long>(v_.size())) {
        v_.resize(static_cast<size_t>(k + 1), -1);
    }
    if (v_[static_cast<size_t>(k)] < 0) ++count_;
    v_[static_cast<size_t>(k)] = idx;
}

int ParamMap::at(int number) const {
    if (!has(number)) throw SnaqError("KeyError: key " + std::to_string(number) + " not found in parameter map");
    return v_[static_cast<size_t>(static_cast<long long>(number) - base_)];
}

// ------------------------------------------------------------------------
// parameters.jl
// ------------------------------------------------------------------------

static bool skip_edge_param(const Network& net, EdgeId e) {
    NodeId child = net.getchild(e);
    if (net.N(child).leaf) return true;
    if (net.N(child).hybrid) {
        auto ch = net.getchildren(child);
        if (ch.size() == 1 && net.N(ch[0]).leaf) return true;
    }
    return false;
}

int countoptimizationparams(const Network& net) {
    int n = static_cast<int>(net.hybrid.size());
    for (EdgeId e : net.edge)
        if (!skip_edge_param(net, e)) ++n;
    return n;
}

OptInfo gatheroptimizationinfo(Network& net, bool change_numbers) {
    OptInfo info;
    if (change_numbers) {
        int uq = net.numedges;
        for (NodeId h : net.hybrid) net.N(h).number = uq++;
        for (EdgeId e : net.edge) net.E(e).number = uq++;
        for (NodeId n : net.node) net.N(n).number = uq++;
    }
    std::vector<ParamObj> objs;
    std::vector<int> numbers;
    for (NodeId h : net.hybrid) {
        objs.push_back({true, h});
        numbers.push_back(net.N(h).number);
    }
    for (EdgeId e : net.edge) {
        objs.push_back({false, e});
        numbers.push_back(net.E(e).number);
    }
    std::vector<size_t> order(objs.size());
    std::iota(order.begin(), order.end(), size_t{0});
    std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) { return numbers[a] < numbers[b]; });
    int param_idx = 0;
    for (size_t k : order) {
        const ParamObj& o = objs[k];
        if (!o.is_hybrid_node && skip_edge_param(net, o.id)) continue;
        int num = numbers[k];
        if (info.param_map.has(num)) throw SnaqError("Duplicate object number #" + std::to_string(num) + ".");
        info.param_map.set(num, param_idx);
        info.idx_obj.push_back(o);
        ++param_idx;
    }
    info.narg = param_idx;
    info.params.resize(static_cast<size_t>(info.narg));
    info.LB.resize(static_cast<size_t>(info.narg));
    info.UB.resize(static_cast<size_t>(info.narg));
    info.init_steps.resize(static_cast<size_t>(info.narg));
    for (int j = 0; j < info.narg; ++j) {
        const ParamObj& o = info.idx_obj[static_cast<size_t>(j)];
        size_t sj = static_cast<size_t>(j);
        if (o.is_hybrid_node) {
            info.params[sj] = net.E(net.getparentedgeminor(o.id)).gamma;
            info.LB[sj] = 0.0;
            info.UB[sj] = 1.0;
            info.init_steps[sj] = 0.1;
        } else {
            info.params[sj] = net.E(o.id).length;
            info.LB[sj] = 0.0;
            info.UB[sj] = 25.0;
            info.init_steps[sj] = 1.0;
        }
    }
    return info;
}

std::vector<double> gatherparams(const Network& net, const ParamMap& pm) {
    std::vector<double> params(pm.size(), 0.0);
    for (NodeId h : net.hybrid)
        if (pm.has(net.N(h).number))
            params[static_cast<size_t>(pm.at(net.N(h).number))] = net.E(net.getparentedgeminor(h)).gamma;
    for (EdgeId e : net.edge)
        if (pm.has(net.E(e).number)) params[static_cast<size_t>(pm.at(net.E(e).number))] = net.E(e).length;
    return params;
}

void setX(Network& net, const double* X, const std::vector<ParamObj>& idx_obj) {
    for (size_t j = 0; j < idx_obj.size(); ++j) {
        const ParamObj& o = idx_obj[j];
        if (o.is_hybrid_node) {
            EdgeId major = net.getparentedge(o.id);
            EdgeId minor = net.getparentedgeminor(o.id);
            net.E(major).gamma = 1.0 - X[j];
            net.E(minor).gamma = X[j];
        } else {
            net.E(o.id).length = X[j];
        }
    }
}

void fixnegativeedges(Network& net, bool shouldwarn) {
    for (EdgeId e : net.edge) {
        if (net.N(net.getchild(e)).leaf) continue;
        Edge& ed = net.E(e);
        if (ed.hybrid && ed.length < 0.0) {
            if (shouldwarn)
                warn("Found a hybrid edge with unspecified length; setting its length to 0.0 before proceeding. To "
                     "avoid this, set the edge's length manually.");
            shouldwarn = false;
            ed.length = 0.0;
        }
        if (ed.length < 0)
            throw SnaqError(
                "Found an edge with unspecified length that was neither a hybrid edge nor a leaf edge. This is not "
                "allowed.\nThis is likely due to the edge not having a specified length in the newick string used to "
                "read the network.");
    }
}

// ------------------------------------------------------------------------
// ancestor paths (ancestor_paths.jl)
// ------------------------------------------------------------------------

namespace {

using EdgeVec = std::vector<EdgeId>;

std::vector<NodeId> ancestorchain(const Network& net, NodeId node) {
    std::vector<NodeId> chain{node};
    while (true) {
        NodeId p = NONE;
        for (EdgeId e : net.N(node).edge)
            if (net.ischildof(node, e)) {
                p = net.getparent(e);
                break;
            }
        if (p == NONE) break;
        node = p;
        chain.push_back(node);
    }
    return chain;
}

std::vector<NodeId> treelikeancestorchain(const Network& net, NodeId node) {
    std::vector<NodeId> chain{node};
    while (true) {
        int npa = 0;
        NodeId p = NONE;
        for (EdgeId e : net.N(node).edge)
            if (net.ischildof(node, e)) {
                if (npa == 0) p = net.getparent(e);
                ++npa;
            }
        if (npa != 1) break;
        if (net.N(p).hybrid) break;
        node = p;
        chain.push_back(node);
    }
    return chain;
}

void pathedges(const Network& net, EdgeVec& dest, const std::vector<NodeId>& cx, const std::vector<NodeId>& cy) {
    size_t ix = cx.size(), iy = cy.size();  // 1-based semantics
    while (ix > 1 && iy > 1 && cx[ix - 2] == cy[iy - 2]) {
        --ix;
        --iy;
    }
    dest.clear();
    for (size_t k = 0; k + 1 < ix; ++k) dest.push_back(net.getparentedge(cx[k]));
    for (size_t k = 0; k + 1 < iy; ++k) dest.push_back(net.getparentedge(cy[k]));
}

bool meetingpathedges(const Network& net, EdgeVec& dest, const std::vector<NodeId>& cx,
                      const std::vector<NodeId>& cy) {
    for (size_t ix = 0; ix < cx.size(); ++ix) {
        NodeId nx = cx[ix];
        for (size_t iy = 0; iy < cy.size(); ++iy) {
            if (cy[iy] != nx) continue;
            dest.clear();
            for (size_t k = 0; k < ix; ++k) dest.push_back(net.getparentedge(cx[k]));
            for (size_t k = 0; k < iy; ++k) dest.push_back(net.getparentedge(cy[k]));
            return true;
        }
    }
    return false;
}

bool contains(const EdgeVec& v, EdgeId e) { return std::find(v.begin(), v.end(), e) != v.end(); }

bool aredisjointedges(const EdgeVec& A, const EdgeVec& B) {
    const EdgeVec& s = A.size() <= B.size() ? A : B;
    const EdgeVec& l = A.size() <= B.size() ? B : A;
    for (EdgeId e : s)
        if (contains(l, e)) return false;
    return true;
}

std::vector<int> common_params(const Network& net, const EdgeVec& A, const EdgeVec& B, const ParamMap& pm) {
    std::vector<int> out;
    for (EdgeId e : A)
        if (contains(B, e)) out.push_back(pm.at(net.E(e).number));
    return out;
}

std::vector<int> unionedgeparams(const Network& net, const EdgeVec& A, const EdgeVec& B, const ParamMap& pm) {
    std::vector<int> out;
    for (EdgeId e : A) {
        int num = net.E(e).number;
        if (!pm.has(num)) continue;
        int p = pm.at(num);
        if (std::find(out.begin(), out.end(), p) == out.end()) out.push_back(p);
    }
    for (EdgeId e : B) {
        if (contains(A, e)) continue;
        int num = net.E(e).number;
        if (!pm.has(num)) continue;
        int p = pm.at(num);
        if (std::find(out.begin(), out.end(), p) == out.end()) out.push_back(p);
    }
    return out;
}

struct PathScratch {
    EdgeVec ab, cd, ac, bd;
};

// treelikequartetdata: returns false when the quartet is not tree-like.
bool treelikequartetdata(const Network& net, const std::vector<NodeId>& ca, const std::vector<NodeId>& cb,
                         const std::vector<NodeId>& cc, const std::vector<NodeId>& cdch, const ParamMap& pm,
                         bool hashybrids, PathScratch& sc, RecursiveCFEquation& eqn, std::vector<int>* relevant) {
    if (hashybrids) {
        if (!meetingpathedges(net, sc.ab, ca, cb)) return false;
        if (!meetingpathedges(net, sc.cd, cc, cdch)) return false;
        if (!meetingpathedges(net, sc.ac, ca, cc)) return false;
        if (!meetingpathedges(net, sc.bd, cb, cdch)) return false;
    } else {
        pathedges(net, sc.ab, ca, cb);
        pathedges(net, sc.cd, cc, cdch);
        pathedges(net, sc.ac, ca, cc);
        pathedges(net, sc.bd, cb, cdch);
    }
    eqn.can_coalesce_here = true;
    eqn.division_H = -1;
    eqn.divisions.clear();
    if (aredisjointedges(sc.ab, sc.cd)) {
        eqn.coal_edges = common_params(net, sc.ac, sc.bd, pm);
        eqn.which_coal = 1;
        if (relevant) *relevant = unionedgeparams(net, sc.ac, sc.bd, pm);
    } else if (aredisjointedges(sc.ac, sc.bd)) {
        eqn.coal_edges = common_params(net, sc.ab, sc.cd, pm);
        eqn.which_coal = 2;
        if (relevant) *relevant = unionedgeparams(net, sc.ab, sc.cd, pm);
    } else {
        eqn.coal_edges = common_params(net, sc.ab, sc.cd, pm);
        eqn.which_coal = 3;
        if (relevant) *relevant = unionedgeparams(net, sc.ab, sc.cd, pm);
    }
    return true;
}

// trytreelikequartet(net, taxa, param_map): network version, by leaf names.
bool trytreelikequartet_names(const Network& net, const std::array<std::string, 4>& taxa, const ParamMap& pm,
                              RecursiveCFEquation& eqn) {
    bool hashybrids = net.numhybrids != 0;
    std::array<std::vector<NodeId>, 4> ch;
    for (int k = 0; k < 4; ++k) {
        NodeId l = NONE;
        for (NodeId x : net.leaf)
            if (net.N(x).name == taxa[static_cast<size_t>(k)]) {
                l = x;
                break;
            }
        if (l == NONE)
            throw SnaqError("trytreelikequartet: taxon " + taxa[static_cast<size_t>(k)] +
                            " is not a leaf of the network.");
        ch[static_cast<size_t>(k)] = hashybrids ? treelikeancestorchain(net, l) : ancestorchain(net, l);
    }
    PathScratch sc;
    return treelikequartetdata(net, ch[0], ch[1], ch[2], ch[3], pm, hashybrids, sc, eqn, nullptr);
}

// ------------------------------------------------------------------------
// hybrid_quartets.jl
// ------------------------------------------------------------------------

Network inducedquartetnetwork(const TreeQuartetContext& ctx, const std::array<int, 4>& taxa) {
    const Network& net = *ctx.net;
    std::unordered_set<NodeId> keep;
    std::vector<NodeId> stack, kept;
    for (int t : taxa) {
        NodeId l = ctx.leaf_of_taxon[static_cast<size_t>(t)];
        if (keep.count(l)) continue;
        keep.insert(l);
        stack.push_back(l);
        kept.push_back(l);
    }
    while (!stack.empty()) {
        NodeId n = stack.back();
        stack.pop_back();
        for (NodeId p : net.getparents(n)) {
            if (keep.count(p)) continue;
            keep.insert(p);
            stack.push_back(p);
            kept.push_back(p);
        }
    }
    std::sort(kept.begin(), kept.end(), [&](NodeId a, NodeId b) {
        return ctx.nodepos[static_cast<size_t>(a)] < ctx.nodepos[static_cast<size_t>(b)];
    });
    Network sub;
    sub.npool.reserve(kept.size());
    std::unordered_map<NodeId, NodeId> nmap;
    for (NodeId n : kept) {
        const Node& src = net.N(n);
        NodeId id = sub.new_node(src.number, src.leaf, src.hybrid);
        sub.N(id).name = src.name;
        nmap[n] = id;
        sub.node.push_back(id);
    }
    std::vector<EdgeId> edges;
    std::unordered_set<EdgeId> seen;
    for (NodeId n : kept)
        for (EdgeId e : net.N(n).edge) {
            if (seen.count(e)) continue;
            if (!(keep.count(net.E(e).nd[0]) && keep.count(net.E(e).nd[1]))) continue;
            seen.insert(e);
            edges.push_back(e);
        }
    std::sort(edges.begin(), edges.end(), [&](EdgeId a, EdgeId b) {
        return ctx.edgepos[static_cast<size_t>(a)] < ctx.edgepos[static_cast<size_t>(b)];
    });
    sub.epool.reserve(edges.size());
    for (EdgeId e : edges) {
        const Edge& src = net.E(e);
        NodeId c1 = nmap.at(src.nd[0]), c2 = nmap.at(src.nd[1]);
        Edge ne;
        ne.number = src.number;
        ne.length = src.length;
        ne.hybrid = src.hybrid;
        ne.gamma = src.hybrid ? src.gamma : 1.0;
        ne.nd[0] = c1;
        ne.nd[1] = c2;
        ne.nnode = 2;
        ne.ischild1 = src.ischild1;
        ne.ismajor = src.ismajor;
        ne.containroot = src.containroot;
        sub.epool.push_back(ne);
        EdgeId id = static_cast<EdgeId>(sub.epool.size() - 1);
        sub.N(c1).edge.push_back(id);
        sub.N(c2).edge.push_back(id);
        sub.edge.push_back(id);
    }
    for (NodeId id : sub.node) {
        if (sub.N(id).leaf) sub.leaf.push_back(id);
        if (sub.N(id).hybrid) sub.hybrid.push_back(id);
    }
    NodeId root = net.getroot();
    if (!keep.count(root)) throw SnaqError("inducedquartetnetwork: the network root is not an ancestor of the quartet.");
    sub.rooti = sub.index_of_node(nmap.at(root));
    sub.numtaxa = static_cast<int>(sub.leaf.size());
    sub.numnodes = static_cast<int>(sub.node.size());
    sub.numedges = static_cast<int>(sub.edge.size());
    sub.numhybrids = static_cast<int>(sub.hybrid.size());
    sub.isrooted = net.isrooted;
    return sub;
}

NodeId getlowesthybridrecur(const Network& net, NodeId node) {
    if (net.N(node).leaf) return NONE;
    for (NodeId child : net.getchildren(node)) {
        NodeId v = getlowesthybridrecur(net, child);
        if (v != NONE) return v;
    }
    return net.N(node).hybrid ? node : NONE;
}

NodeId getlowesthybrid(const Network& net) {
    if (net.numhybrids == 1) return net.hybrid[0];
    return getlowesthybridrecur(net, net.hybrid[0]);
}

int nleavesbelowlowesthybrid(const Network& net, NodeId H) {
    std::vector<NodeId> queue = net.getchildren(H);
    int found = 0;
    while (!queue.empty()) {
        NodeId curr = queue.back();
        queue.pop_back();
        if (net.N(curr).leaf)
            ++found;
        else
            for (NodeId c : net.getchildren(curr)) queue.push_back(c);
    }
    return found;
}

std::vector<EdgeId> getinternaledgesbelowlowesthybrid(const Network& net, NodeId H) {
    std::vector<EdgeId> out;
    std::vector<NodeId> c = net.getchildren(H);
    while (c.size() == 1) {
        out.push_back(net.getparentedge(c[0]));
        c = net.getchildren(c[0]);
    }
    return out;
}

std::vector<NodeId> getleavesbelowlowesthybrid(const Network& net, NodeId H) {
    std::vector<NodeId> queue{H}, leaves;
    while (!queue.empty()) {
        NodeId curr = queue.back();
        queue.pop_back();
        if (net.N(curr).leaf) leaves.push_back(curr);
        for (NodeId c : net.getchildren(curr)) queue.push_back(c);
    }
    return leaves;
}

NodeId find_hybrid(const Network& net, const std::string& name, int number) {
    for (NodeId h : net.hybrid)
        if (net.N(h).name == name && net.N(h).number == number) return h;
    throw SnaqError("hybrid node not found in network copy");
}

NodeId find_leaf(const Network& net, const std::string& name) {
    for (NodeId l : net.leaf)
        if (net.N(l).name == name) return l;
    throw SnaqError("leaf " + name + " not found in network copy");
}

const DeleteOpts kKeepRootNoFuse{/*nofuse*/ true, /*simplify*/ false, /*unroot*/ false, /*multgammas*/ false,
                                 /*keeporiginalroot*/ true};

// Remove edge `rm` (a parent edge of hybrid H) and its references in `div`,
// as in the n_below_H == 1 case of get4taxaquartetequations.
void dropparentedge(Network& div, EdgeId rm, EdgeId keepe) {
    NodeId other_parent = div.getparent(rm);
    for (int k = 0; k < div.E(rm).nnode; ++k) {
        NodeId nd = div.E(rm).nd[k];
        int i = div.N(nd).edge.find(rm);
        // node.edge = [e for e in node.edge if e != rm]
        while (i >= 0) {
            div.N(nd).edge.erase_at(static_cast<size_t>(i));
            i = div.N(nd).edge.find(rm);
        }
    }
    div.deleteEdge(rm);
    div.removeHybrid(div.getchild(keepe));
    if (div.getchildren(other_parent).empty()) {
        div.N(other_parent).leaf = true;
        div.leaf.push_back(other_parent);
        deleteleaf_node(div, other_parent, kKeepRootNoFuse);
        div.N(other_parent).leaf = false;
    }
    div.E(keepe).hybrid = false;
    div.E(keepe).ismajor = true;
    div.N(div.getchild(keepe)).hybrid = false;
}

RecursiveCFEquation get4taxaquartetequations(Network& net, const std::array<std::string, 4>& taxa,
                                             const ParamMap& pm) {
    if (net.numhybrids == 0) {
        RecursiveCFEquation eqn;
        trytreelikequartet_names(net, taxa, pm, eqn);
        return eqn;
    }
    // Special case: we may still have hybrids, but 3+ leaves share a parent
    for (size_t k = 0; k < 2; ++k) {
        NodeId L = net.leaf.at(k);
        NodeId p_leaf = net.getparent_node(L);
        auto ch = net.getchildren(p_leaf);
        int nleaves = 0;
        for (NodeId c : ch) nleaves += net.N(c).leaf ? 1 : 0;
        if (ch.size() > 2 && nleaves > 2) {
            RecursiveCFEquation eqn;
            eqn.can_coalesce_here = true;
            eqn.which_coal = 1;
            eqn.division_H = -1;
            return eqn;
        }
    }
    NodeId lowest_H = getlowesthybrid(net);
    int n_below_H = nleavesbelowlowesthybrid(net, lowest_H);
    const std::string Hname = net.N(lowest_H).name;
    const int Hnumber = net.N(lowest_H).number;
    if (n_below_H == 1 &&
        net.getparent(net.getparentedge(lowest_H)) != net.getparent(net.getparentedgeminor(lowest_H))) {
        Network div_major = net.deepcopynetwork();
        {
            NodeId H = find_hybrid(div_major, Hname, Hnumber);
            EdgeId E_minor = div_major.getparentedgeminor(H);
            EdgeId E_major = div_major.getparentedge(H);
            dropparentedge(div_major, E_minor, E_major);
        }
        Network div_minor = net.deepcopynetwork();
        {
            NodeId H = find_hybrid(div_minor, Hname, Hnumber);
            EdgeId E_minor = div_minor.getparentedgeminor(H);
            EdgeId E_major = div_minor.getparentedge(H);
            dropparentedge(div_minor, E_major, E_minor);
        }
        RecursiveCFEquation out;
        out.can_coalesce_here = false;
        out.which_coal = 0;
        out.division_H = pm.at(Hnumber);
        out.divisions.reserve(2);
        out.divisions.push_back(get4taxaquartetequations(div_minor, taxa, pm));
        out.divisions.push_back(get4taxaquartetequations(div_major, taxa, pm));
        return out;
    } else if (n_below_H == 2) {
        std::vector<EdgeId> int_edges = getinternaledgesbelowlowesthybrid(net, lowest_H);
        std::vector<NodeId> leaves_below_H = getleavesbelowlowesthybrid(net, lowest_H);
        std::vector<std::string> below_names;
        for (NodeId L : leaves_below_H) below_names.push_back(net.N(L).name);
        std::array<std::string, 2> leaf_names{below_names[0], below_names[1]};
        std::sort(leaf_names.begin(), leaf_names.end());

        // Both taxa take the minor edge (div1) / the major edge (div2)
        auto both_take = [&](bool take_minor) {
            Network div = net.deepcopynetwork();
            NodeId H = find_hybrid(div, Hname, Hnumber);
            EdgeId E_minor = div.getparentedgeminor(H);
            EdgeId E_major = div.getparentedge(H);
            EdgeId take = take_minor ? E_minor : E_major;
            EdgeId drop = take_minor ? E_major : E_minor;
            for (const std::string& lname : below_names) {
                NodeId divL = find_leaf(div, lname);
                NodeId l = addleaf(div, div.getchild(take), "__" + lname, 0.0);
                deleteleaf_node(div, divL, kKeepRootNoFuse);
                div.N(l).name = lname;
            }
            deletehybridedge(div, drop, true, false, false, true, true);
            return div;
        };
        // Lower taxon takes `first`, higher takes the other (div3: minor/major; div4: major/minor)
        auto split_take = [&](bool lower_takes_minor) {
            Network div = net.deepcopynetwork();
            NodeId H = find_hybrid(div, Hname, Hnumber);
            EdgeId E_minor = div.getparentedgeminor(H);
            EdgeId E_major = div.getparentedge(H);
            EdgeId e_low = lower_takes_minor ? E_minor : E_major;
            EdgeId e_high = lower_takes_minor ? E_major : E_minor;
            NodeId nl = addleaf(div, div.getparent(e_low), "__" + leaf_names[0], 0.0);
            deleteleaf_node(div, find_leaf(div, leaf_names[0]), kKeepRootNoFuse);
            div.N(nl).name = leaf_names[0];
            nl = addleaf(div, div.getparent(e_high), "__" + leaf_names[1], 0.0);
            deleteleaf_node(div, find_leaf(div, leaf_names[1]), kKeepRootNoFuse);
            div.N(nl).name = leaf_names[1];
            return div;
        };
        Network div1 = both_take(true);
        Network div2 = both_take(false);
        Network div3 = split_take(true);
        Network div4 = split_take(false);

        int which_quartet;
        if (leaf_names[0] == taxa[0])
            which_quartet = (leaf_names[1] == taxa[1]) ? 1 : (leaf_names[1] == taxa[2] ? 2 : 3);
        else if (leaf_names[0] == taxa[1])
            which_quartet = (leaf_names[1] == taxa[2]) ? 3 : 2;
        else
            which_quartet = 1;

        RecursiveCFEquation out;
        out.divisions.reserve(4);
        out.divisions.push_back(get4taxaquartetequations(div1, taxa, pm));
        out.divisions.push_back(get4taxaquartetequations(div2, taxa, pm));
        out.divisions.push_back(get4taxaquartetequations(div3, taxa, pm));
        out.divisions.push_back(get4taxaquartetequations(div4, taxa, pm));
        out.can_coalesce_here = !int_edges.empty();
        for (EdgeId e : int_edges) out.coal_edges.push_back(pm.at(net.E(e).number));
        out.which_coal = which_quartet;
        out.division_H = pm.at(Hnumber);
        return out;
    } else {
        // 3 or 4 leaves below this hybrid, so it has no effect on eCFs
        deletehybridedge(net, net.getparentedgeminor(lowest_H), false, true, false, true, false);
        return get4taxaquartetequations(net, taxa, pm);
    }
}

QuartetData reticulatequartetequations(const TreeQuartetContext& ctx, const std::array<int, 4>& taxa,
                                       const ParamMap& pm) {
    Network net = inducedquartetnetwork(ctx, taxa);
    auto bcc = biconnectedcomponents(net, true);
    auto entry = biconnectedcomponent_entrynodes(net, bcc, true);
    std::vector<long> entryindex;
    for (NodeId n : entry) {
        auto it = std::find(net.vec_node.begin(), net.vec_node.end(), n);
        entryindex.push_back(static_cast<long>(it - net.vec_node.begin()));
    }
    auto exitnodes = biconnectedcomponent_exitnodes(net, bcc, false);
    std::vector<size_t> bloborder(bcc.size());
    std::iota(bloborder.begin(), bloborder.end(), size_t{0});
    std::stable_sort(bloborder.begin(), bloborder.end(),
                     [&](size_t a, size_t b) { return entryindex[a] < entryindex[b]; });
    auto isexternal = [&](size_t ib) {
        if (exitnodes[ib].size() != 1) return false;
        auto ch = net.getchildren(exitnodes[ib][0]);
        return ch.size() == 1 && net.N(ch[0]).leaf;
    };
    for (auto it = bloborder.rbegin(); it != bloborder.rend(); ++it) {
        size_t ib = *it;
        if (!isexternal(ib)) continue;
        for (EdgeId he : bcc[ib]) {
            if (net.E(he).ismajor) continue;
            if (net.index_of_edge(he) < 0) continue;
            deletehybridedge(net, he, false, true, false, true, false);
        }
    }
    std::array<std::string, 4> names;
    for (int k = 0; k < 4; ++k) names[static_cast<size_t>(k)] = ctx.taxa[static_cast<size_t>(taxa[static_cast<size_t>(k)])];
    auto q = std::make_shared<QuartetDataImpl>();
    q->eqn = get4taxaquartetequations(net, names, pm);
    for (EdgeId e : net.edge)
        if (pm.has(net.E(e).number)) q->relevant_params.push_back(pm.at(net.E(e).number));
    for (NodeId h : net.hybrid)
        if (pm.has(net.N(h).number)) q->relevant_params.push_back(pm.at(net.N(h).number));
    q->q_taxa = taxa;
    return q;
}

}  // namespace

// ------------------------------------------------------------------------
// contexts and batches
// ------------------------------------------------------------------------

TreeQuartetContext treequartetcontext(const Network& net) {
    TreeQuartetContext ctx;
    ctx.net = &net;
    ctx.hashybrids = net.numhybrids != 0;
    std::vector<std::pair<std::string, NodeId>> lv;
    for (NodeId l : net.leaf) lv.emplace_back(net.N(l).name, l);
    std::sort(lv.begin(), lv.end());
    for (auto& p : lv) {
        ctx.taxa.push_back(p.first);
        ctx.leaf_of_taxon.push_back(p.second);
        ctx.chains.push_back(ctx.hashybrids ? treelikeancestorchain(net, p.second) : ancestorchain(net, p.second));
    }
    if (ctx.hashybrids) {
        ctx.nodepos.assign(net.npool.size(), -1);
        ctx.edgepos.assign(net.epool.size(), -1);
        for (size_t j = 0; j < net.node.size(); ++j) ctx.nodepos[static_cast<size_t>(net.node[j])] = static_cast<int>(j);
        for (size_t j = 0; j < net.edge.size(); ++j) ctx.edgepos[static_cast<size_t>(net.edge[j])] = static_cast<int>(j);
    }
    return ctx;
}

static thread_local PathScratch tl_scratch;

QuartetData findquartetequations4taxa(const TreeQuartetContext& ctx, const std::array<int, 4>& taxa,
                                      const ParamMap& pm) {
    const auto& ch = ctx.chains;
    auto q = std::make_shared<QuartetDataImpl>();
    if (treelikequartetdata(*ctx.net, ch[static_cast<size_t>(taxa[0])], ch[static_cast<size_t>(taxa[1])],
                            ch[static_cast<size_t>(taxa[2])], ch[static_cast<size_t>(taxa[3])], pm, ctx.hashybrids,
                            tl_scratch, q->eqn, &q->relevant_params)) {
        q->q_taxa = taxa;
        return q;
    }
    return reticulatequartetequations(ctx, taxa, pm);
}

static void check_gammas(const Network& net) {
    for (EdgeId e : net.edge)
        if (net.E(e).hybrid && !(net.E(e).gamma <= 1.0 && net.E(e).gamma >= 0.0))
            throw SnaqError("net has gammas that are not in [0, 1]");
    for (NodeId h : net.hybrid) {
        double s = net.E(net.getparentedge(h)).gamma + net.E(net.getparentedgeminor(h)).gamma;
        if (!(std::abs(s - 1.0) <= 1.4901161193847656e-8 * std::max(std::abs(s), 1.0)))
            throw SnaqError("net has hybrid with gammas that do not sum to 1");
    }
}

EquationSet findquartetequations(Network& net, const std::vector<long long>& ranks, int nthreads) {
    fixnegativeedges(net);
    check_gammas(net);
    EquationSet out;
    out.info = gatheroptimizationinfo(net, true);
    TreeQuartetContext ctx = treequartetcontext(net);
    out.taxa = ctx.taxa;
    const int ntax = static_cast<int>(ctx.taxa.size());
    out.eqns.resize(ranks.size());
    parallel_for(ranks.size(), nthreads, [&](size_t lo, size_t hi, int) {
        for (size_t j = lo; j < hi; ++j)
            out.eqns[j] = findquartetequations4taxa(ctx, unrank4taxa(ntax, ranks[j]), out.info.param_map);
    });
    return out;
}

EquationSet findquartetequations(Network& net, int nthreads) {
    long long nq = nchoose4(net.numtaxa);
    std::vector<long long> ranks(static_cast<size_t>(nq));
    std::iota(ranks.begin(), ranks.end(), 0LL);
    return findquartetequations(net, ranks, nthreads);
}

void rebuildquartetequations(Network& net, const std::vector<QuartetData>& old_eqns,
                             std::vector<QuartetData>& new_eqns, int nthreads) {
    OptInfo info = gatheroptimizationinfo(net, true);
    TreeQuartetContext ctx = treequartetcontext(net);
    new_eqns.resize(old_eqns.size());
    parallel_for(old_eqns.size(), nthreads, [&](size_t lo, size_t hi, int) {
        for (size_t j = lo; j < hi; ++j)
            new_eqns[j] = findquartetequations4taxa(ctx, old_eqns[j]->q_taxa, info.param_map);
    });
}

// ------------------------------------------------------------------------
// in-place updates
// ------------------------------------------------------------------------

static std::vector<int> paramsbelowurNNI1(const Network& net, NodeId u, const ParamMap& pm) {
    std::vector<int> objs_below;
    auto in_objs = [&](int x) { return std::find(objs_below.begin(), objs_below.end(), x) != objs_below.end(); };
    if (!net.getparents(u).empty()) {
        objs_below.push_back(pm.at(net.E(net.getparentedge(u)).number));
    } else {
        for (EdgeId e : net.N(u).edge)
            if (pm.has(net.E(e).number)) objs_below.push_back(pm.at(net.E(e).number));
    }
    std::vector<EdgeId> queue;
    for (EdgeId e : net.N(u).edge)
        if (net.getparent(e) == u) queue.push_back(e);
    long long n_iter = 0;
    while (!queue.empty()) {
        if (++n_iter >= 10000000) throw SnaqError("Stuck in infinite loop in paramsbelowurNNI1.");
        EdgeId next_edge = queue.back();
        queue.pop_back();
        int num = net.E(next_edge).number;
        if (pm.has(num)) {
            // QUIRK: compares an edge *number* against parameter indices.
            if (in_objs(num))
                continue;
            else
                objs_below.push_back(pm.at(num));
        }
        NodeId child = net.getchild(next_edge);
        // QUIRK: same number-vs-index comparison as above.
        if (net.N(child).hybrid && !in_objs(net.N(child).number)) objs_below.push_back(pm.at(net.N(child).number));
        for (EdgeId e : net.N(child).edge)
            if (e != next_edge && net.getchild(e) != child) queue.push_back(e);
    }
    return objs_below;
}

static bool recurfxnhasparams(const RecursiveCFEquation& q, const std::vector<int>& idxs) {
    for (int p : idxs) {
        if (std::find(q.coal_edges.begin(), q.coal_edges.end(), p) != q.coal_edges.end()) return true;
        if (q.division_H == p) return true;
    }
    for (const auto& d : q.divisions)
        if (recurfxnhasparams(d, idxs)) return true;
    return false;
}

void applyrNNI1update(Network& Nprime, const std::vector<QuartetData>& old_q, std::vector<QuartetData>& new_q,
                      const ParamMap& pm, NodeId u, int nthreads) {
    std::vector<int> relevant = paramsbelowurNNI1(Nprime, u, pm);
    TreeQuartetContext ctx = treequartetcontext(Nprime);
    new_q.resize(old_q.size());
    parallel_for(old_q.size(), nthreads, [&](size_t lo, size_t hi, int) {
        for (size_t j = lo; j < hi; ++j) {
            const auto& o = old_q[j];
            if (!o->eqn.divisions.empty() || recurfxnhasparams(o->eqn, relevant))
                new_q[j] = findquartetequations4taxa(ctx, o->q_taxa, pm);
            else
                new_q[j] = o;
        }
    });
}

void applyrNNI2update(Network& Nprime, const std::vector<QuartetData>& old_q, std::vector<QuartetData>& new_q,
                      const ParamMap& pm, NodeId s, NodeId t, NodeId u, NodeId v, int nthreads) {
    auto edge_to = [&](NodeId a, NodeId b) {
        for (EdgeId e : Nprime.N(a).edge)
            if (Nprime.E(e).has_node(b)) return e;
        throw SnaqError("applyrNNI2update: nodes are not adjacent");
    };
    std::vector<int> relevant{pm.at(Nprime.E(edge_to(u, t)).number), pm.at(Nprime.E(edge_to(s, v)).number)};
    TreeQuartetContext ctx = treequartetcontext(Nprime);
    new_q.resize(old_q.size());
    parallel_for(old_q.size(), nthreads, [&](size_t lo, size_t hi, int) {
        for (size_t j = lo; j < hi; ++j) {
            const auto& o = old_q[j];
            bool contains_param = false;
            for (int p : relevant)
                if (std::find(o->relevant_params.begin(), o->relevant_params.end(), p) != o->relevant_params.end())
                    contains_param = true;
            new_q[j] = contains_param ? findquartetequations4taxa(ctx, o->q_taxa, pm) : o;
        }
    });
}

}  // namespace csnaq
