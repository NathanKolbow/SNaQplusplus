#include "newick.hpp"

#include <cctype>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>

#include "network_manip.hpp"

namespace csnaq {

namespace {

constexpr int END = -1;  // "missing" character (end of input)

bool is_space(int c) { return c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\f' || c == '\v'; }

bool isvalidnewick_namesymbol(int c) {
    if (c == END) return false;
    if (std::isalnum(c) || c == '_' || c == '#') return true;
    return !is_space(c) && c != '(' && c != ')' && c != '[' && c != ']' && c != ':' && c != ';' && c != ',';
}

class Stream {
public:
    explicit Stream(const std::string& s) : s_(s) {}
    bool eof() const { return pos_ >= s_.size(); }
    // peekskip: skip whitespace, return next char without consuming it
    int peekskip() {
        while (pos_ < s_.size() && is_space(static_cast<unsigned char>(s_[pos_]))) ++pos_;
        return pos_ < s_.size() ? static_cast<unsigned char>(s_[pos_]) : END;
    }
    // readskip!: read chars until a non-space one, return the last read
    int readskip() {
        int c = END;
        while (pos_ < s_.size()) {
            c = static_cast<unsigned char>(s_[pos_++]);
            if (!is_space(c)) break;
        }
        return c;
    }
    // read(s, Char)
    int read() { return pos_ < s_.size() ? static_cast<unsigned char>(s_[pos_++]) : END; }
    std::string rest() const { return s_.substr(std::min(pos_, s_.size())); }

private:
    const std::string& s_;
    size_t pos_ = 0;
};

struct Parser {
    Stream s;
    Network& net;
    int numLeft = 1;
    std::vector<std::string> hybrids;

    Parser(const std::string& str, Network& n) : s(str), net(n) {}

    [[noreturn]] void fail(const std::string& msg) { throw SnaqError(msg); }

    int advance() {
        int c = s.readskip();
        if (c == END)
            fail("Tree ends prematurely while reading subtree after left parenthesis " + std::to_string(numLeft - 1) +
                 ".");
        return c;
    }

    void readnexus_comment(int c) {
        while (c == '[') {
            if (s.read() != '[') fail("I was supposed to read '['");
            if (s.read() != '&') fail("read '[' but not followed by &");
            int d;
            do {
                d = s.read();
            } while (d != END && d != ']');
            if (d == END) fail("comment without ] to end it");
            c = s.peekskip();
        }
    }

    struct NameResult {
        int num;
        std::string name;
        bool pound;
    };

    NameResult readnewick_nodename(int c) {
        int pound = 0;
        std::string name;
        while (isvalidnewick_namesymbol(c)) {
            s.readskip();  // advance past c
            if (c == '#') {
                ++pound;
                c = s.readskip();
                if (!(c != END && std::isalpha(c)))
                    fail("Expected name after # but received " + std::string(1, static_cast<char>(c)) +
                         " in left parenthesis " + std::to_string(numLeft - 1) + ". remaining: " + s.rest() + ".");
                if (c != 'L' && c != 'H' && c != 'R')
                    warn("Expected H, R or LGT after # but received " + std::string(1, static_cast<char>(c)) +
                         " in left parenthesis " + std::to_string(numLeft - 1) + ".");
                name = std::string(1, static_cast<char>(c)) + name;
            } else {
                name.push_back(static_cast<char>(c));
            }
            c = s.peekskip();
        }
        if (pound > 1) fail("strange node name with " + std::to_string(pound) + " # signs: " + name);
        readnexus_comment(c);
        return {static_cast<int>(net.names.size()) + 1, name, pound == 1};
    }

    double readnewick_float(int c) {
        auto isfloatchar = [](int ch) {
            return (ch != END && std::isdigit(ch)) || ch == '.' || ch == 'e' || ch == '-' || ch == 'E' || ch == '+';
        };
        if (!isfloatchar(c)) fail("Expected float digit after ':' but found " + std::string(1, static_cast<char>(c)));
        std::string num;
        while (isfloatchar(c)) {
            num.push_back(static_cast<char>(s.read()));
            c = s.peekskip();
        }
        try {
            size_t idx = 0;
            double v = std::stod(num, &idx);
            if (idx != num.size()) throw std::invalid_argument(num);
            return v;
        } catch (...) {
            fail("problem with number read " + num + ", not a float number");
        }
    }

    double getfloat() {
        int c = s.peekskip();
        if (c == '[') {
            readnexus_comment(c);
            c = s.peekskip();
        }
        if (c == ':' || c == ')' || c == ',') return -1.0;
        double val = readnewick_float(c);
        if (val < 0.0) {
            std::cerr << "┌ Error: expecting non-negative value but read '-', left parenthesis " << (numLeft - 1)
                      << ". will set to 0.\n";
            val = 0.0;
        }
        readnexus_comment(s.peekskip());
        return val;
    }

    void edgedata(EdgeId e) {
        s.read();  // the first ':'
        net.E(e).length = getfloat();
        if (s.peekskip() == ':') {
            s.readskip();
            net.E(e).y = getfloat();
        }
        if (s.peekskip() == ':') {
            s.readskip();
            net.E(e).gamma = getfloat();
        }
        Edge& ed = net.E(e);
        if (ed.gamma != -1.0 && !ed.hybrid && ed.gamma != 1.0) {
            warn("γ read for edge " + std::to_string(ed.number) + " but it is not hybrid, so γ=" +
                 julia_float(ed.gamma) + " ignored");
            ed.gamma = 1.0;
        }
    }

    void synchronizepartnersdata(EdgeId e, NodeId n) {
        std::vector<EdgeId> partners;
        for (EdgeId e2 : net.N(n).edge)
            if (net.E(e2).hybrid && e2 != e && n == net.getchild(e2)) partners.push_back(e2);
        if (partners.empty()) return;
        if (partners.size() > 1)
            fail("γ value parsed but hybrid edge has " + std::to_string(partners.size()) +
                 " partners (should be 0 or 1)");
        Edge& ed = net.E(e);
        Edge& partner = net.E(partners[0]);
        if (ed.gamma == -1.0 && partner.gamma == -1.0) return;
        if (ed.gamma == -1.0) {
            if (partner.gamma < 0.0) {
                warn("partners: edge with no γ, partner with γ<0. will turn to 1 and 0");
                partner.gamma = 0.0;
                ed.gamma = 1.0;
            } else if (partner.gamma > 1.0) {
                warn("partners: edge with no γ, partner with γ>1. will turn to 0 and 1");
                partner.gamma = 1.0;
                ed.gamma = 0.0;
            } else {
                ed.gamma = 1.0 - partner.gamma;
            }
        } else if (partner.gamma == -1.0) {
            if (ed.gamma < 0.0) {
                warn("partners: edge with no γ, partner with γ<0. will turn to 1 and 0");
                ed.gamma = 0.0;
                partner.gamma = 1.0;
            } else if (ed.gamma > 1.0) {
                warn("partners: edge with no γ, partner with γ>1. will turn to 0 and 1");
                ed.gamma = 1.0;
                partner.gamma = 0.0;
            } else {
                partner.gamma = 1.0 - ed.gamma;
            }
        } else {
            double gammasum = ed.gamma + partner.gamma;
            // isapprox(gammasum, 1.0): rtol = sqrt(eps)
            if (!(std::abs(gammasum - 1.0) <= 1.4901161193847656e-8 * std::max(std::abs(gammasum), 1.0))) {
                ed.gamma = ed.gamma / gammasum;
                partner.gamma = 1.0 - ed.gamma;
            }
        }
        if (!(std::abs(ed.gamma - 0.5) <= 1.4901161193847656e-8 * std::max(std::abs(ed.gamma), 0.5))) {
            bool emajor = ed.gamma > 0.5;
            if (ed.ismajor != emajor) {
                ed.ismajor = emajor;
                partner.ismajor = !emajor;
            }
        }
    }

    // parsenewick_hybridnode!: returns (edge, retained node)
    std::pair<EdgeId, NodeId> hybridnode(NodeId n, NodeId parent, const std::string& name) {
        net.N(n).hybrid = true;
        bool seen = false;
        for (const auto& h : hybrids)
            if (h == name) {
                seen = true;
                break;
            }
        EdgeId e = net.new_edge(net.numedges + 1);
        if (net.N(n).leaf) net.E(e).ismajor = false;
        net.E(e).hybrid = true;
        net.E(e).gamma = -1.0;
        if (seen) {
            NodeId other = NONE;
            for (NodeId no : net.node)
                if (net.N(no).name == name) {
                    other = no;
                    break;
                }
            if (other == NONE) fail("hybrid name " + name + " was supposed to be in the network, but not found");
            if (!net.N(n).leaf && !net.N(other).leaf) {
                fail("both hybrid nodes are internal nodes: successors of the hybrid node must only be included "
                     "in the node list of a single occurrence of the hybrid node.");
            } else if (net.N(n).leaf) {
                net.pushEdge(e);
                net.setNode(e, other, parent);
                net.setEdge(other, e);
                net.setEdge(parent, e);
                n = other;
            } else {
                if (net.N(other).edge.size() != 1)
                    fail("strange: node " + std::to_string(net.N(other).number) +
                         " is a leaf hybrid node. should have only 1 edge");
                EdgeId otheredge = net.N(other).edge[0];
                net.removeNode(other, otheredge);
                net.deleteNode(other);
                net.setNode(otheredge, n);
                net.setEdge(n, otheredge);
                net.setNode(e, n, parent);
                net.setEdge(n, e);
                net.setEdge(parent, e);
                net.pushNode(n);
                net.pushEdge(e);
                net.N(n).number = net.N(other).number;
                net.N(n).name = net.N(other).name;
            }
        } else {
            net.N(n).hybrid = true;
            net.names.push_back(name);
            net.N(n).name = name;
            hybrids.push_back(name);
            net.pushNode(n);
            net.pushEdge(e);
            net.setNode(e, n, parent);
            net.setEdge(n, e);
            net.setEdge(parent, e);
        }
        net.E(e).containroot = !net.E(e).hybrid;
        return {e, n};
    }

    EdgeId treenode(NodeId n, NodeId parent) {
        net.pushNode(n);
        EdgeId e = net.new_edge(net.numedges + 1);
        net.pushEdge(e);
        net.setNode(e, n, parent);
        net.setEdge(n, e);
        net.setEdge(parent, e);
        return e;
    }

    NodeId remainingsubtree() {
        numLeft += 1;
        NodeId n = net.new_node(-numLeft, false);
        s.readskip();  // the '('
        bool keepon = true;
        while (keepon) {
            subtree(n);
            int c = advance();
            if (c == ')')
                keepon = false;
            else if (c != ',')
                fail("Expected right parenthesis after left parenthesis " + std::to_string(numLeft - 1) +
                     " but read " + std::string(1, static_cast<char>(c)) + ". The remainder of line is " + s.rest() +
                     ".");
        }
        return n;
    }

    void subtree(NodeId parent) {
        int c = s.peekskip();
        bool hasname = false;
        NameResult nr;
        NodeId n;
        if (c == '(') {
            n = remainingsubtree();
            c = s.peekskip();
            nr = readnewick_nodename(c);
            if (!nr.name.empty()) {
                hasname = true;
                net.N(n).number = nr.num;
            }
        } else {
            hasname = true;
            nr = readnewick_nodename(c);
            if (nr.name.empty())
                fail("Expected digit, alphanum or # at the start of taxon name, but received " +
                     (c == END ? std::string("end of input") : std::string(1, static_cast<char>(c))) +
                     ". remaining: " + s.rest() + ".");
            n = net.new_node(nr.num, true);
        }
        EdgeId e;
        if (nr.pound) {
            auto r = hybridnode(n, parent, nr.name);
            e = r.first;
            n = r.second;
        } else {
            if (hasname) {
                net.names.push_back(nr.name);
                net.N(n).name = nr.name;
            }
            e = treenode(n, parent);
        }
        c = s.peekskip();
        net.E(e).length = -1.0;
        if (c == ':') edgedata(e);
        if (net.E(e).hybrid) synchronizepartnersdata(e, n);
    }

    void checkNumHybEdges() {
        if (net.isTree()) return;
        if (net.hybrid.empty()) fail("net.hybrid should not be empty for this network");
        for (NodeId n : net.hybrid) {
            int hyb = 0;
            for (EdgeId e : net.N(n).edge) hyb += net.E(e).hybrid ? 1 : 0;
            if (hyb == 1) {
                fail("hybrid node " + std::to_string(net.N(n).number) + " named " + net.N(n).name +
                     " has only one hybrid edge attached");
            } else if (hyb == 0) {
                size_t ne = net.N(n).edge.size();
                if (ne == 0)
                    fail("strange hybrid node attached to 0 edges");
                else if (ne == 1) {
                    if (!net.N(n).leaf) fail("hybrid node has only one tree edge attached and it is not a leaf");
                } else {
                    warn("hybrid node " + std::to_string(net.N(n).number) +
                         " not connected to any hybrid edges. Now transformed to tree node");
                    net.N(n).hybrid = false;
                }
            } else {
                int nhybparents = 0;
                for (EdgeId e : net.N(n).edge)
                    if (n == net.getchild(e)) {
                        if (net.E(e).hybrid)
                            ++nhybparents;
                        else
                            std::cerr << "┌ Error: node " << net.N(n).number
                                      << " has parent tree edge: wrong ischild1 for this edge?\n";
                    }
                if (nhybparents < 2) fail("hybrid node with fewer than 2 hybrid parents");
                if (nhybparents > 2) fail("hybrid node with more than 2 hybrid parents: we don't allow such polytomies");
            }
        }
    }

    void parse() {
        // readuntil(s, ";"): the description must contain a ';'
        // (checked by the caller)
        int c = s.peekskip();
        if (c != '(') fail("Expected beginning of tree with ( but received something else instead, rest is " + s.rest());
        numLeft += 1;
        NodeId n = net.new_node(-numLeft, false);
        c = s.readskip();
        while (c != ';') {
            subtree(n);
            c = s.readskip();
            if (s.eof() && c != ';') {
                fail("Tree ended while reading in subtree beginning with left parenthesis number " +
                     std::to_string(numLeft - 1) + ".");
            } else if (c == ',') {
                continue;
            } else if (c == ')') {
                c = s.peekskip();
                NameResult nr = readnewick_nodename(c);
                if (!nr.name.empty()) net.N(n).name = nr.name;
                c = s.peekskip();
                if (c == ':') {
                    while (c != ';' && c != END) c = s.readskip();
                }
                if (c == END) fail("file does not end in ;");
            }
        }
        if (net.N(n).edge.size() == 1) {
            EdgeId edge = net.N(n).edge[0];
            NodeId child = net.getOtherNode(edge, n);
            net.removeEdge(child, edge);
            net.rooti = net.index_of_node(child);
            net.deleteEdge(edge);
        } else {
            net.pushNode(n);
            net.rooti = net.index_of_node(n);
        }
        // storeHybrids!
        std::vector<NodeId> hy;
        for (NodeId x : net.node)
            if (net.N(x).hybrid) hy.push_back(x);
        if (!hy.empty()) {
            net.hybrid = hy;
            net.numhybrids = static_cast<int>(hy.size());
        }
        checkNumHybEdges();
        directedges(net, true);
        net.isrooted = true;
    }
};

std::string lstrip(const std::string& s) {
    size_t i = 0;
    while (i < s.size() && is_space(static_cast<unsigned char>(s[i]))) ++i;
    return s.substr(i);
}

}  // namespace

Network readnewick_string(const std::string& input) {
    std::string str = lstrip(input);
    if (str.empty() || str[0] != '(') throw SnaqError("Newick description must start with '('");
    // Julia: readuntil(s, ";", keep=true) must end with ';'
    if (str.find(';') == std::string::npos) throw SnaqError("file does not end in ;");
    Network net;
    Parser p(str, net);
    p.parse();
    return net;
}

std::vector<Network> readmultinewick_file(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw SnaqError("Could not find or open " + path + " file");
    std::vector<Network> out;
    std::string line;
    size_t lineno = 0;
    while (std::getline(in, line)) {
        ++lineno;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        std::string t = lstrip(line);
        if (t.empty()) continue;
        try {
            out.push_back(readnewick_string(t));
        } catch (const SnaqError& e) {
            throw SnaqError("error reading line " + std::to_string(lineno) + " of " + path + ": " + e.what());
        }
    }
    return out;
}

std::vector<Network> read_networks_arg(const std::string& arg) {
    std::string t = lstrip(arg);
    if (!t.empty() && t[0] == '(') return {readnewick_string(t)};
    return readmultinewick_file(arg);
}

// ------------------------------------------------------------------------
// writing
// ------------------------------------------------------------------------

namespace {

std::string fmt(double x, bool roundBL, int digits) { return julia_float(roundBL ? julia_round(x, digits) : x); }

void writesubtree_node(std::string& s, const Network& net, NodeId n, EdgeId parent, bool di, bool roundBL, int digits,
                       bool internallabel) {
    const Node& nd = net.N(n);
    if (!nd.leaf && (parent == NONE || net.E(parent).ismajor)) {
        s += "(";
        bool firstchild = true;
        for (EdgeId e : nd.edge) {
            if (e == parent) continue;
            if (parent == NONE && n == net.getchild(e)) continue;
            if (net.E(e).hybrid && net.getchild(e) == n) continue;
            if (!firstchild) s += ",";
            firstchild = false;
            writesubtree_node(s, net, net.getOtherNode(e, n), e, di, roundBL, digits, internallabel);
        }
        s += ")";
    }
    if (parent != NONE && net.E(parent).hybrid) {
        s += "#";
        s += nd.name;
        if (nd.name.empty() && !net.E(parent).ismajor)
            warn("hybrid node " + std::to_string(nd.number) + " has no name");
    } else if (internallabel || nd.leaf) {
        s += nd.name;
    }
    bool printBL = false;
    if (parent != NONE && net.E(parent).length != -1.0) {
        s += ":";
        s += fmt(net.E(parent).length, roundBL, digits);
        printBL = true;
    }
    if (parent != NONE && net.E(parent).hybrid && !di) {
        double g = net.E(parent).gamma;
        if (g != -1.0) {
            s += (printBL ? "::" : ":::");
            s += fmt(g, roundBL, digits);
        }
    }
    if (parent == NONE) s += ";";
}

}  // namespace

std::string writenewick(Network& net, bool round, int digits, bool di, bool internallabel) {
    int rootsaved = net.rooti;
    bool changeroot = false;
    try {
        directedges(net);
    } catch (const RootMismatch& err) {
        std::cout << err.what() << "\nCannot write topology with current root." << std::endl;
        changeroot = true;
    }
    while (changeroot) {
        for (EdgeId e : net.edge) {
            int i = net.index_of_node(net.getparent(e));
            net.rooti = i;
            try {
                directedges(net);
                std::cout << "Setting root at node " << net.N(net.node[static_cast<size_t>(i)]).number
                          << " (net.rooti = " << (i + 1) << ")\n\n";
                changeroot = false;
                break;
            } catch (const RootMismatch&) {
            }
        }
        if (changeroot) {
            net.rooti = rootsaved;
            throw RootMismatch("Could not find admissible root. Cannot write topology.");
        }
    }
    if (net.N(net.getroot()).leaf)
        warn("Root is placed at a leaf node, so the parenthetical format will look strange.\nUse rootatnode! or "
             "rootonedge! to change the root position");
    std::string s;
    NodeId rootnode = net.getroot();
    if (net.numnodes > 1) {
        s += "(";
        size_t degree = net.N(rootnode).edge.size();
        for (EdgeId e : net.N(rootnode).edge) {
            writesubtree_node(s, net, net.getOtherNode(e, rootnode), e, di, round, digits, internallabel);
            --degree;
            if (degree != 0) s += ",";
        }
        s += ")";
    }
    if (internallabel || net.numnodes == 1) s += net.N(rootnode).name;
    s += ";";
    return s;
}

}  // namespace csnaq
