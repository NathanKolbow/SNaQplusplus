#include "test_support.hpp"

#include <algorithm>
#include <atomic>
#include <filesystem>
#include <fstream>
#include <set>
#include <sstream>

#ifndef CSNAQ_TEST_DATA_DIR
#error "CSNAQ_TEST_DATA_DIR must be defined"
#endif

namespace csnaq::test {

namespace fs = std::filesystem;

std::string data(const std::string& file) { return std::string(CSNAQ_TEST_DATA_DIR) + "/" + file; }

std::string tempdir() {
    static std::atomic<int> counter{0};
    fs::path p = fs::temp_directory_path() /
                 ("csnaq_test_" + std::to_string(static_cast<long long>(wall_time() * 1e6)) + "_" +
                  std::to_string(counter++));
    fs::create_directories(p);
    return p.string();
}

Network nw(const std::string& newick) { return readnewick_string(newick); }

std::vector<Network> read_trees(const std::string& fixture) { return readmultinewick_file(data(fixture)); }

std::vector<std::string> read_lines(const std::string& fixture) {
    std::ifstream in(data(fixture));
    std::vector<std::string> out;
    std::string l;
    while (std::getline(in, l)) out.push_back(l);
    return out;
}

CFMatrix read_cf(const std::string& fixture) {
    std::ifstream in(data(fixture));
    size_t n = 0;
    in >> n;
    CFMatrix m;
    m.v.resize(3 * n);
    for (size_t i = 0; i < 3 * n; ++i) in >> m.v[i];
    return m;
}

NodeId node_named(const Network& net, const std::string& name) {
    for (NodeId n : net.node)
        if (net.N(n).name == name) return n;
    throw SnaqError("no node named " + name);
}

std::vector<NodeId> get_nodes(const Network& net, std::initializer_list<const char*> names) {
    std::vector<NodeId> out;
    for (const char* n : names) out.push_back(node_named(net, n));
    return out;
}

bool testnewickequality(const Network& net, const std::string& correct) {
    Network c = net.clone();
    for (EdgeId e : c.edge) c.E(e).length = -1;
    return writenewick(c) == correct;
}

// ------------------------------------------------------------------------
// hardwired clusters (PhyloNetworks compareNetworks.jl)
// ------------------------------------------------------------------------

namespace {

using Row = std::vector<int>;  // [edge number, cluster indicators..., (type)]

int find_taxon(const std::vector<std::string>& S, const std::string& name) {
    auto it = std::find(S.begin(), S.end(), name);
    if (it == S.end()) throw SnaqError("leaf " + name + " not in taxon list");
    return static_cast<int>(it - S.begin());
}

void traverseTree2Matrix(const Network& T, NodeId node, EdgeId edge, std::vector<Row>& M, size_t& ie,
                         const std::vector<std::string>& S) {
    NodeId child = T.getOtherNode(edge, node);
    size_t indedge = ie;
    M[indedge][0] = T.E(edge).number;
    ++ie;
    for (EdgeId e : T.N(child).edge) {
        if (e == edge) continue;
        NodeId gc = T.getOtherNode(e, child);
        if (T.N(gc).leaf) {
            M[indedge][static_cast<size_t>(find_taxon(S, T.N(gc).name)) + 1] = 1;
        } else {
            size_t inde = ie;
            traverseTree2Matrix(T, child, e, M, ie, S);
            for (size_t k = 1; k < M[indedge].size(); ++k) M[indedge][k] |= M[inde][k];
        }
    }
}

std::vector<Row> tree2Matrix(const Network& T, const std::vector<std::string>& S, bool rooted) {
    size_t ne = T.edge.size() - static_cast<size_t>(T.numtaxa);
    NodeId root = T.getroot();
    if (T.N(root).leaf) ++ne;
    std::vector<Row> M(ne, Row(S.size() + 1, 0));
    size_t ie = 0;
    for (EdgeId e : T.N(root).edge)
        if (!T.N(T.getOtherNode(e, root)).leaf) traverseTree2Matrix(T, root, e, M, ie, S);
    if (!rooted && T.N(root).edge.size() < 3 && !M.empty()) M.erase(M.begin());
    return M;
}

int hardwiredclusters_rec(const Network& net, NodeId node, EdgeId edge, std::vector<Row>& M, size_t& ie,
                          const std::vector<std::string>& S) {
    NodeId child = net.getOtherNode(edge, node);
    if (net.N(child).leaf) {
        find_taxon(S, net.N(child).name);
        return -1;
    }
    EdgeId partner = NONE;
    if (net.E(edge).hybrid) {
        if (net.getchild(edge) != child) throw SnaqError("inconsistency during network traversal");
        bool partnerVisited = true;
        int indpartner = -1;
        for (EdgeId e : net.N(child).edge)
            if (net.E(e).hybrid && e != edge && net.getchild(e) == child) {
                partner = e;
                for (size_t r = 0; r < ie; ++r)
                    if (M[r][0] == net.E(e).number) {
                        indpartner = static_cast<int>(r);
                        break;
                    }
                if (indpartner < 0) partnerVisited = false;
                break;
            }
        if (partner == NONE) throw SnaqError("partner hybrid edge not found");
        if (partnerVisited) return indpartner;
    }
    size_t indedge = ie;
    M[indedge][0] = net.E(edge).number;
    M[indedge].back() = net.E(edge).hybrid ? 11 : 10;
    ++ie;
    for (EdgeId e : net.N(child).edge) {
        if (e == edge || (net.E(edge).hybrid && e == partner)) continue;
        NodeId gc = net.getOtherNode(e, child);
        if (net.N(gc).leaf) {
            M[indedge][static_cast<size_t>(find_taxon(S, net.N(gc).name)) + 1] = 1;
        } else {
            int inde = hardwiredclusters_rec(net, child, e, M, ie, S);
            for (size_t k = 1; k + 1 < M[indedge].size(); ++k) M[indedge][k] |= M[static_cast<size_t>(inde)][k];
        }
    }
    return static_cast<int>(indedge);
}

std::vector<Row> hardwiredclusters(const Network& net, const std::vector<std::string>& S) {
    size_t ne = net.edge.size() - static_cast<size_t>(net.numtaxa) - net.hybrid.size();
    NodeId root = net.getroot();
    if (net.N(root).leaf) ++ne;
    std::vector<Row> M(ne, Row(S.size() + 2, 0));
    size_t ie = 0;
    for (EdgeId e : net.N(root).edge) hardwiredclusters_rec(net, root, e, M, ie, S);
    return M;
}

int hwcd_rooted_or_trees(const Network& net1, const Network& net2, bool rooted) {
    bool bothtrees = net1.numhybrids == 0 && net2.numhybrids == 0;
    std::vector<std::string> taxa = net1.tiplabels();
    std::sort(taxa.begin(), taxa.end());
    auto t2 = net2.tiplabels();
    std::set<std::string> s2(t2.begin(), t2.end());
    for (const auto& t : taxa)
        if (!s2.count(t)) throw SnaqError("net1 and net2 do not share the same taxon set.");
    std::vector<Row> M1, M2;
    if (bothtrees) {
        M1 = tree2Matrix(net1, taxa, rooted);
        M2 = tree2Matrix(net2, taxa, rooted);
    } else {
        M1 = hardwiredclusters(net1, taxa);
        M2 = hardwiredclusters(net2, taxa);
    }
    int dis = 0;
    std::vector<size_t> n2ci(M2.size());
    for (size_t k = 0; k < M2.size(); ++k) n2ci[k] = k;
    for (const Row& r1 : M1) {
        Row a(r1.begin() + 1, r1.end());
        Row m1(a.size());
        for (size_t k = 0; k < a.size(); ++k) m1[k] = 1 - a[k];
        bool found = false;
        for (size_t j = n2ci.size(); j-- > 0;) {
            const Row& r2 = M2[n2ci[j]];
            Row b(r2.begin() + 1, r2.end());
            if (a == b || (!rooted && m1 == b)) {
                found = true;
                n2ci.erase(n2ci.begin() + static_cast<long>(j));
                break;
            }
        }
        if (!found) ++dis;
    }
    return dis + dis + static_cast<int>(M2.size()) - static_cast<int>(M1.size());
}

}  // namespace

int hardwiredclusterdistance(const Network& n1, const Network& n2, bool rooted) {
    bool bothtrees = n1.numhybrids == 0 && n2.numhybrids == 0;
    if (rooted || bothtrees) return hwcd_rooted_or_trees(n1, n2, rooted);
    Network net1 = n1.clone(), net2 = n2.clone();
    removedegree2nodes(net1);
    removedegree2nodes(net2);
    auto roots = [](Network& net) {
        std::vector<int> r;
        for (NodeId n : net.node)
            if (!net.N(n).leaf) r.push_back(net.N(n).number);
        for (size_t i = r.size(); i-- > 0;) {
            try {
                rootatnode_number(net, r[i]);
            } catch (const RootMismatch&) {
                r.erase(r.begin() + static_cast<long>(i));
            }
        }
        return r;
    };
    std::vector<int> r1 = roots(net1), r2 = roots(net2);
    int best = std::numeric_limits<int>::max();
    for (int a : r1) {
        rootatnode_number(net1, a);
        for (int b : r2) {
            rootatnode_number(net2, b);
            best = std::min(best, hwcd_rooted_or_trees(net1, net2, true));
        }
    }
    return best;
}

bool shrinkedge(Network& net, EdgeId e) {
    if (net.E(e).hybrid) throw SnaqError("cannot shrink hybrid edge");
    NodeId cn = net.getchild(e);
    if (net.N(cn).hybrid) throw SnaqError("cannot shrink tree edge: its child node is a hybrid");
    NodeId pn = net.getparent(e);
    if (net.isexternal(e)) throw SnaqError("won't shrink edge: it is incident to a leaf");
    net.removeEdge(pn, e);
    net.E(e).nnode = 0;
    for (EdgeId ee : net.N(cn).edge) {
        if (ee == e) continue;
        Edge& ed = net.E(ee);
        for (int k = 0; k < ed.nnode; ++k)
            if (ed.nd[k] == cn) {
                ed.nd[k] = pn;
                break;
            }
        net.N(pn).edge.push_back(ee);
    }
    bool anyh = false;
    for (EdgeId x : net.N(pn).edge) anyh |= net.E(x).hybrid;
    net.N(pn).booln1 = anyh;
    net.N(cn).edge.clear();
    net.deleteEdge(e);
    net.deleteNode(cn);
    if (net.N(pn).hybrid) {
        int nc = 0;
        for (EdgeId x : net.N(pn).edge) nc += (!net.E(x).hybrid || net.getchild(x) != pn) ? 1 : 0;
        return nc > 1;
    }
    return false;
}

// ------------------------------------------------------------------------
// tcgidentifiable (identifiability_properties.jl)
// ------------------------------------------------------------------------

namespace {

std::vector<NodeId> getleavesbelow(const Network& net, std::vector<NodeId> nodes) {
    std::vector<NodeId> leaves, visited;
    while (!nodes.empty()) {
        NodeId curr = nodes.front();
        nodes.erase(nodes.begin());
        if (std::find(visited.begin(), visited.end(), curr) != visited.end()) continue;
        visited.push_back(curr);
        if (net.N(curr).leaf) {
            leaves.push_back(curr);
            continue;
        }
        for (NodeId c : net.getchildren(curr)) nodes.push_back(c);
    }
    return leaves;
}

std::vector<int> getdiamondns(const Network& net, const std::vector<NodeId>& cycle) {
    NodeId n3 = NONE, n1 = NONE;
    for (NodeId n : cycle)
        if (net.N(n).hybrid) {
            n3 = n;
            break;
        }
    for (NodeId n : cycle) {
        if (n == n3) continue;
        bool adj = false;
        for (EdgeId e : net.N(n).edge) adj |= net.E(e).has_node(n3);
        if (!adj) {
            n1 = n;
            break;
        }
    }
    std::vector<NodeId> rest;
    for (NodeId n : cycle)
        if (n != n3 && n != n1) rest.push_back(n);
    std::vector<NodeId> roots{rest.at(0), n1, rest.at(1), n3};
    std::vector<int> ns(4, 0);
    for (size_t k = 0; k < 4; ++k) {
        std::vector<NodeId> visited = roots;
        std::vector<NodeId> Q{roots[k]};
        while (!Q.empty()) {
            NodeId curr = Q.front();
            Q.erase(Q.begin());
            visited.push_back(curr);
            if (net.N(curr).leaf) {
                ++ns[k];
            } else {
                for (EdgeId e : net.N(curr).edge) {
                    NodeId a = net.getOtherNode(e, curr);
                    if (std::find(visited.begin(), visited.end(), a) != visited.end()) continue;
                    Q.push_back(a);
                }
            }
        }
    }
    return ns;
}

std::pair<bool, int> isCk(const Network& N0) {
    Network N = N0.deepcopynetwork();
    if (!isgallednetwork(N)) return {false, 0};
    if (!std::get<2>(istreechild(N))) return {false, 0};
    if (shrink2cycles(N)) return {false, 0};
    int k = N.numtaxa;
    auto bcc = biconnectedcomponents(N);
    for (const auto& blob : bcc) {
        if (blob.size() <= 1) continue;
        std::vector<NodeId> lnodes;
        for (EdgeId e : blob) {
            NodeId c = N.getchild(e);
            if (N.N(c).hybrid && std::find(lnodes.begin(), lnodes.end(), c) == lnodes.end()) lnodes.push_back(c);
        }
        std::vector<NodeId> Y = getleavesbelow(N, lnodes);
        for (size_t l = 0; l < Y.size(); ++l) {
            Network bloblet = N.deepcopynetwork();
            for (size_t i = 0; i < Y.size(); ++i) {
                if (i == l) continue;
                deleteleaf_name(bloblet, N.N(Y[i]).name, DeleteOpts{false, true, false, false, false});
            }
            auto bb = biconnectedcomponents(bloblet);
            std::vector<size_t> sizes;
            for (const auto& c : bb)
                if (c.size() > 1) sizes.push_back(c.size());
            if (sizes.size() != 1) return {false, 0};
            k = std::min(k, static_cast<int>(sizes[0]));
        }
    }
    return {true, k};
}

}  // namespace

bool tcgidentifiable(const Network& N0) {
    if (N0.numhybrids == 0) return true;
    if (N0.isrooted) return false;
    Network N = N0.deepcopynetwork();
    if (shrink2cycles(N) || shrink3cycles(N)) return false;
    if (getnetworklevel(N) == 1) {
        auto bcc = biconnectedcomponents(N);
        for (const auto& blob : bcc) {
            if (blob.size() <= 1) continue;
            std::vector<NodeId> cyc;
            for (EdgeId e : blob)
                for (NodeId n : {N.getchild(e), N.getparent(e)})
                    if (std::find(cyc.begin(), cyc.end(), n) == cyc.end()) cyc.push_back(n);
            if (cyc.size() != 4) continue;
            auto ns = getdiamondns(N, cyc);
            if (!(ns[0] >= 2 || ns[2] >= 2 || (ns[1] >= 2 && ns[3] >= 2))) return false;
        }
    } else {
        auto [ok, k] = isCk(N);
        if (!ok) return false;
        if (k < 4) return false;
        if (k > 4) return true;
        for (NodeId ln : N.hybrid)
            if (getleavesbelow(N, {ln}).size() < 2) return false;
    }
    return true;
}

// ------------------------------------------------------------------------

Network random_tree(int n, Rng& rng) {
    // random binary tree by sequential random joins, internal nodes named
    std::vector<std::string> parts;
    for (int j = 1; j <= n; ++j) parts.push_back("t" + std::to_string(j));
    while (parts.size() > 1) {
        size_t a = rng.index(parts.size());
        std::string x = parts[a];
        parts.erase(parts.begin() + static_cast<long>(a));
        size_t b = rng.index(parts.size());
        std::string y = parts[b];
        parts.erase(parts.begin() + static_cast<long>(b));
        double l1 = 0.1 + rng.rand(), l2 = 0.1 + rng.rand();
        parts.push_back("(" + x + ":" + julia_float(l1) + "," + y + ":" + julia_float(l2) + ")i" +
                        std::to_string(rng.randint(0, 9999)));
    }
    Network t = nw(parts[0] + ";");
    // make internal names unique
    int k = 0;
    for (NodeId nd : t.node)
        if (!t.N(nd).leaf) t.N(nd).name = "i" + std::to_string(++k);
    return t;
}

SnaqOptions quiet_options() {
    SnaqOptions o;
    o.filename = "";
    o.threads = 4;
    o.seed = 42;
    o.seed_given = true;
    return o;
}

}  // namespace csnaq::test
