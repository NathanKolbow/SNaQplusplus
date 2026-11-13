#include "quartets.hpp"

#include <algorithm>
#include <cmath>
#include <set>
#include <thread>
#include <unordered_map>

namespace csnaq {

std::array<int, 4> unrank4taxa(int ntaxa, long long rank) {
    if (rank < 0 || rank >= nchoose4(ntaxa))
        throw SnaqError("the index for the quartet " + std::to_string(rank + 1) +
                        " needs to be less than choose(n,4)=" + std::to_string(nchoose4(ntaxa)));
    std::array<int, 4> out{};
    long long r = rank;
    int upper = ntaxa - 1;
    for (int k = 4; k >= 1; --k) {
        // largest t <= upper with C(t, k) <= r (binary search; C(k-1, k) = 0 <= r)
        int lo = k - 1, hi = upper;
        while (lo < hi) {
            int mid = lo + (hi - lo + 1) / 2;
            if (binomial(mid, k) <= r)
                lo = mid;
            else
                hi = mid - 1;
        }
        int t = lo;
        out[static_cast<size_t>(k - 1)] = t;
        r -= binomial(t, k);
        upper = t - 1;
    }
    return out;
}

long long rank4taxa(const std::array<int, 4>& s) {
    return binomial(s[0], 1) + binomial(s[1], 2) + binomial(s[2], 3) + binomial(s[3], 4);
}

long long quartetsamplesize(long long n, double p) {
    long long a = static_cast<long long>(std::ceil(static_cast<double>(n) * p));
    return std::max(a, std::min(n, 10LL));
}

std::vector<long long> sampleqindices(long long n, double p, Rng& rng) {
    std::vector<long long> s = rng.sample_without_replacement(n, quartetsamplesize(n, p));
    std::sort(s.begin(), s.end());
    return s;
}

std::vector<long long> sampleqindices(long long n, double p, const std::vector<char>& informative, Rng& rng) {
    std::vector<long long> valid;
    for (size_t i = 0; i < informative.size(); ++i)
        if (informative[i]) valid.push_back(static_cast<long long>(i));
    long long ninform = static_cast<long long>(valid.size());
    long long nsamp = std::min(ninform, static_cast<long long>(std::ceil(static_cast<double>(n) * p)));
    if (ninform == nsamp) return valid;
    std::vector<long long> pick = rng.sample_without_replacement(ninform, nsamp);
    std::vector<long long> out;
    out.reserve(pick.size());
    for (long long k : pick) out.push_back(valid[static_cast<size_t>(k)]);
    std::sort(out.begin(), out.end());
    return out;
}

bool isquartetinformative(const double* row, double atol) {
    for (int i = 0; i < 2; ++i)
        for (int j = i + 1; j < 3; ++j)
            if (std::abs(row[i] - row[j]) > atol) return true;
    return false;
}

std::vector<std::string> sorted_union_taxa(const std::vector<Network>& trees) {
    std::set<std::string> s;
    for (const Network& t : trees)
        for (NodeId l : t.leaf) s.insert(t.N(l).name);
    return {s.begin(), s.end()};
}

// ------------------------------------------------------------------------

long long ObservedCFs::pairindex(int i, int j) const {
    const long long n = ntaxa();
    // 0-based version of SNaQ's pairindex
    return static_cast<long long>(i) * (2 * n - i - 1) / 2 + (j - i - 1);
}

ObservedCFs::ObservedCFs(const std::vector<Network>& trees, std::vector<std::string> taxa, bool lazy_semantics)
    : taxa_(std::move(taxa)), ntrees_(static_cast<int>(trees.size())), lazy_(lazy_semantics) {
    std::sort(taxa_.begin(), taxa_.end());
    const int n = ntaxa();
    if (n < 4) throw SnaqError("LazyQuartetCF needs at least 4 taxa (got " + std::to_string(n) + ").");
    for (const Network& t : trees)
        if (t.numhybrids != 0) throw SnaqError("gene tree found that is a network");
    std::unordered_map<std::string, int> taxonindex;
    for (int i = 0; i < n; ++i) taxonindex[taxa_[static_cast<size_t>(i)]] = i;
    const long long npairs = static_cast<long long>(n) * (n - 1) / 2;
    lcadepth_.assign(static_cast<size_t>(npairs * ntrees_), static_cast<int16_t>(-1));

    for (int t = 0; t < ntrees_; ++t) {
        const Network& tre = trees[static_cast<size_t>(t)];
        std::vector<int> buf(tre.leaf.size() + 1);
        // iterative post-order walk replicating pairlcadepths!
        struct Frame {
            NodeId node;
            int depth;
            int start;
            int n;  // taxa collected so far below node
            std::vector<NodeId> children;
            size_t next;
        };
        std::vector<Frame> stack;
        auto push = [&](NodeId nd, int depth, int start) {
            NodeList ch = tre.getchildren(nd);
            Frame f{nd, depth, start, 0, std::vector<NodeId>(ch.begin(), ch.end()), 0};
            stack.push_back(std::move(f));
        };
        push(tre.getroot(), 0, 0);
        int returned = -1;  // number of taxa written by the frame just popped
        while (!stack.empty()) {
            Frame& f = stack.back();
            if (tre.N(f.node).leaf) {
                auto it = taxonindex.find(tre.N(f.node).name);
                if (it == taxonindex.end())
                    throw SnaqError("gene tree leaf " + tre.N(f.node).name + " is not in the taxon list.");
                buf[static_cast<size_t>(f.start)] = it->second;
                returned = 1;
                stack.pop_back();
                continue;
            }
            if (returned >= 0) {
                // a child just finished, writing `returned` taxa at f.start + f.n
                const int m = returned;
                const int16_t d16 = static_cast<int16_t>(f.depth);
                for (int x = f.start; x < f.start + f.n; ++x) {
                    int bx = buf[static_cast<size_t>(x)];
                    for (int y = f.start + f.n; y < f.start + f.n + m; ++y) {
                        int by = buf[static_cast<size_t>(y)];
                        int a = std::min(bx, by), b = std::max(bx, by);
                        lcadepth_[static_cast<size_t>(pairindex(a, b) * ntrees_ + t)] = d16;
                    }
                }
                f.n += m;
                returned = -1;
            }
            if (f.next < f.children.size()) {
                NodeId ch = f.children[f.next++];
                int depth = f.depth + 1, start = f.start + f.n;
                if (static_cast<size_t>(start) >= buf.size()) buf.resize(static_cast<size_t>(start) + 1);
                push(ch, depth, start);
            } else {
                returned = f.n;
                stack.pop_back();
            }
        }
    }
}

std::array<double, 3> ObservedCFs::compute(long long rank) const {
    std::array<int, 4> q = unrank4taxa(ntaxa(), rank);
    const int a = q[0], b = q[1], c = q[2], d = q[3];
    const int16_t* pab = &lcadepth_[static_cast<size_t>(pairindex(a, b) * ntrees_)];
    const int16_t* pcd = &lcadepth_[static_cast<size_t>(pairindex(c, d) * ntrees_)];
    const int16_t* pac = &lcadepth_[static_cast<size_t>(pairindex(a, c) * ntrees_)];
    const int16_t* pbd = &lcadepth_[static_cast<size_t>(pairindex(b, d) * ntrees_)];
    const int16_t* pad = &lcadepth_[static_cast<size_t>(pairindex(a, d) * ntrees_)];
    const int16_t* pbc = &lcadepth_[static_cast<size_t>(pairindex(b, c) * ntrees_)];
    long long n1 = 0, n2 = 0, n3 = 0;
    for (int t = 0; t < ntrees_; ++t) {
        int16_t dab = pab[t];
        if (dab < 0) continue;
        int16_t dcd = pcd[t];
        if (dcd < 0) continue;
        int16_t dac = pac[t], dbd = pbd[t], dad = pad[t], dbc = pbc[t];
        int m1 = std::max(dab, dcd), m2 = std::max(dac, dbd), m3 = std::max(dad, dbc);
        if (lazy_) {
            // QUIRK: SNaQ's lazy DataCF counts a quartet that is unresolved in a
            // gene tree (tied LCA depths, i.e. a polytomy) for the first tied topology.
            if (m1 >= m2 && m1 >= m3)
                ++n1;
            else if (m2 >= m3)
                ++n2;
            else
                ++n3;
        } else {
            // Non-lazy DataCF (readtrees2CF): unresolved quartets are skipped.
            if (m1 > m2 && m1 > m3)
                ++n1;
            else if (m2 > m1 && m2 > m3)
                ++n2;
            else if (m3 > m1 && m3 > m2)
                ++n3;
        }
    }
    double suma = static_cast<double>(n1 + n2 + n3);
    return {static_cast<double>(n1) / suma, static_cast<double>(n2) / suma, static_cast<double>(n3) / suma};
}

void ObservedCFs::compute_rows(const std::vector<long long>& ranks, double* out, int nthreads) const {
    const size_t nq = ranks.size();
    auto work = [&](size_t lo, size_t hi) {
        for (size_t j = lo; j < hi; ++j) {
            auto cf = compute(ranks[j]);
            out[3 * j] = cf[0];
            out[3 * j + 1] = cf[1];
            out[3 * j + 2] = cf[2];
        }
    };
    if (nthreads <= 1 || nq < 1024) {
        work(0, nq);
        return;
    }
    std::vector<std::thread> th;
    size_t chunk = (nq + static_cast<size_t>(nthreads) - 1) / static_cast<size_t>(nthreads);
    for (int k = 0; k < nthreads; ++k) {
        size_t lo = static_cast<size_t>(k) * chunk, hi = std::min(nq, lo + chunk);
        if (lo >= hi) break;
        th.emplace_back(work, lo, hi);
    }
    for (auto& t : th) t.join();
}

}  // namespace csnaq
