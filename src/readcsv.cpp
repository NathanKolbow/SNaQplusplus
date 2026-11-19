#include "readcsv.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "quartets.hpp"

namespace snaqpp {

namespace {

// Strip leading/trailing whitespace and optional enclosing quotes from s.
std::string strip(const std::string& s) {
    size_t a = 0, b = s.size();
    while (a < b && (s[a] == ' ' || s[a] == '\t' || s[a] == '\r' || s[a] == '\n')) ++a;
    while (b > a && (s[b - 1] == ' ' || s[b - 1] == '\t' || s[b - 1] == '\r' || s[b - 1] == '\n')) --b;
    if (b >= a + 2 && ((s[a] == '"' && s[b - 1] == '"') || (s[a] == '\'' && s[b - 1] == '\''))) {
        ++a;
        --b;
        while (a < b && (s[a] == ' ' || s[a] == '\t' || s[a] == '\r' || s[a] == '\n')) ++a;
        while (b > a && (s[b - 1] == ' ' || s[b - 1] == '\t' || s[b - 1] == '\r' || s[b - 1] == '\n')) --b;
    }
    return s.substr(a, b - a);
}

// Split a line by comma into fields.
std::vector<std::string> splitcsv(const std::string& line) {
    std::vector<std::string> fields;
    std::istringstream ss(line);
    std::string field;
    while (std::getline(ss, field, ',')) fields.push_back(strip(field));
    return fields;
}

double parse_cf_double(const std::string& s, int lineno, int col) {
    try {
        size_t idx = 0;
        double d = std::stod(s, &idx);
        if (idx != s.size()) throw std::invalid_argument(s);
        return d;
    } catch (...) {
        throw SnaqError("cannot parse CF value at line " + std::to_string(lineno) + ", column " +
                        std::to_string(col + 1) + ": \"" + s + "\"");
    }
}

// Case-insensitive equality.
bool ieq(const std::string& a, const std::string& b) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i)
        if (std::tolower(static_cast<unsigned char>(a[i])) != std::tolower(static_cast<unsigned char>(b[i])))
            return false;
    return true;
}

}  // namespace

CsvCF readtablecsv(const std::string& filename) {
    std::ifstream in(filename);
    if (!in) throw SnaqError("cannot open CF table: " + filename);

    // Read header line.
    std::string header;
    while (std::getline(in, header)) {
        std::string h = strip(header);
        if (!h.empty()) break;
    }
    if (header.empty()) throw SnaqError("CF table is empty: " + filename);
    // Strip UTF-8 BOM if present.
    if (header.size() >= 3 && static_cast<unsigned char>(header[0]) == 0xEF &&
        static_cast<unsigned char>(header[1]) == 0xBB &&
        static_cast<unsigned char>(header[2]) == 0xBF) {
        header = header.substr(3);
    }
    std::vector<std::string> hdr = splitcsv(header);

    // Identify column indices by header names (case-insensitive).
    int col_t1 = -1, col_t2 = -1, col_t3 = -1, col_t4 = -1;
    int col_cf12 = -1, col_cf13 = -1, col_cf14 = -1;
    int col_ngenes = -1;
    for (int i = 0; i < static_cast<int>(hdr.size()); ++i) {
        const std::string& h = hdr[static_cast<size_t>(i)];
        if (ieq(h, "taxon1") || ieq(h, "t1") || ieq(h, "tx1") || ieq(h, "tax1"))
            col_t1 = i;
        else if (ieq(h, "taxon2") || ieq(h, "t2") || ieq(h, "tx2") || ieq(h, "tax2"))
            col_t2 = i;
        else if (ieq(h, "taxon3") || ieq(h, "t3") || ieq(h, "tx3") || ieq(h, "tax3"))
            col_t3 = i;
        else if (ieq(h, "taxon4") || ieq(h, "t4") || ieq(h, "tx4") || ieq(h, "tax4"))
            col_t4 = i;
        else if (ieq(h, "CF12_34") || ieq(h, "CF12.34") || ieq(h, "CF1234") || ieq(h, "obsCF12"))
            col_cf12 = i;
        else if (ieq(h, "CF13_24") || ieq(h, "CF13.24") || ieq(h, "CF1324") || ieq(h, "obsCF13"))
            col_cf13 = i;
        else if (ieq(h, "CF14_23") || ieq(h, "CF14.23") || ieq(h, "CF1423") || ieq(h, "obsCF14"))
            col_cf14 = i;
        else if (ieq(h, "ngenes") || ieq(h, "genes") || ieq(h, "numgenes") || ieq(h, "n_genes"))
            col_ngenes = i;
    }
    if (col_t1 < 0 || col_t2 < 0 || col_t3 < 0 || col_t4 < 0)
        throw SnaqError("CF table " + filename +
                        " must have columns taxon1, taxon2, taxon3, taxon4 (found: " + header + ")");
    if (col_cf12 < 0 || col_cf13 < 0 || col_cf14 < 0)
        throw SnaqError("CF table " + filename +
                        " must have columns CF12_34, CF13_24, CF14_23 (found: " + header + ")");

    const int ncols = static_cast<int>(hdr.size());

    // First pass: collect taxa and all rows.
    struct Row {
        std::string t1, t2, t3, t4;
        double cf12, cf13, cf14;
        double ngenes;
    };
    std::vector<Row> rows;
    std::set<std::string> taxon_set;
    std::string line;
    int lineno = 1;
    while (std::getline(in, line)) {
        ++lineno;
        std::string trimmed = strip(line);
        if (trimmed.empty()) continue;
        std::vector<std::string> fields = splitcsv(trimmed);
        if (static_cast<int>(fields.size()) < ncols)
            throw SnaqError("CF table " + filename + " line " + std::to_string(lineno) + ": expected " +
                            std::to_string(ncols) + " columns, got " + std::to_string(fields.size()));
        Row r;
        r.t1 = fields[static_cast<size_t>(col_t1)];
        r.t2 = fields[static_cast<size_t>(col_t2)];
        r.t3 = fields[static_cast<size_t>(col_t3)];
        r.t4 = fields[static_cast<size_t>(col_t4)];
        r.cf12 = parse_cf_double(fields[static_cast<size_t>(col_cf12)], lineno, col_cf12);
        r.cf13 = parse_cf_double(fields[static_cast<size_t>(col_cf13)], lineno, col_cf13);
        r.cf14 = parse_cf_double(fields[static_cast<size_t>(col_cf14)], lineno, col_cf14);
        r.ngenes = (col_ngenes >= 0) ? parse_cf_double(fields[static_cast<size_t>(col_ngenes)], lineno, col_ngenes) : 0;
        taxon_set.insert(r.t1);
        taxon_set.insert(r.t2);
        taxon_set.insert(r.t3);
        taxon_set.insert(r.t4);
        rows.push_back(std::move(r));
    }
    if (rows.empty()) throw SnaqError("CF table " + filename + " has no data rows");

    // Sorted taxa list and index map.
    std::vector<std::string> taxa(taxon_set.begin(), taxon_set.end());
    std::sort(taxa.begin(), taxa.end());
    const int ntaxa = static_cast<int>(taxa.size());
    if (ntaxa < 4) throw SnaqError("CF table " + filename + " has fewer than 4 taxa");
    std::unordered_map<std::string, int> taxon_index;
    for (int i = 0; i < ntaxa; ++i) taxon_index[taxa[static_cast<size_t>(i)]] = i;

    const long long nq = nchoose4(static_cast<long long>(ntaxa));
    CFMatrix cfs;
    cfs.v.assign(static_cast<size_t>(3 * nq), 1.0 / 3.0);  // default: uninformative
    // Track which quartets have been set for weighting purposes.
    std::vector<bool> present(static_cast<size_t>(nq), false);
    std::vector<double> ngenes_vec;
    if (col_ngenes >= 0) ngenes_vec.assign(static_cast<size_t>(nq), 0.0);

    double max_ngenes = 0.0;

    for (const Row& r : rows) {
        auto it1 = taxon_index.find(r.t1);
        auto it2 = taxon_index.find(r.t2);
        auto it3 = taxon_index.find(r.t3);
        auto it4 = taxon_index.find(r.t4);
        // Should always succeed since we built taxon_index from the same data.
        int i1 = it1->second, i2 = it2->second, i3 = it3->second, i4 = it4->second;
        if (i1 == i2 || i1 == i3 || i1 == i4 || i2 == i3 || i2 == i4 || i3 == i4)
            throw SnaqError("CF table " + filename + ": quartet (" + r.t1 + ", " + r.t2 + ", " + r.t3 + ", " + r.t4 +
                            ") has duplicate taxa");

        // Sort the four taxa indices and figure out how the input CF columns
        // map to the sorted quartet's CF columns.
        // The CSV's CF12_34 is the split {t1,t2}|{t3,t4}. After sorting the
        // four indices into s[0]<s[1]<s[2]<s[3], the three resolutions in
        // SNaQ's order are:
        //   CF(s0,s1 | s2,s3)  = "12|34" of the sorted quartet
        //   CF(s0,s2 | s1,s3)  = "13|24" of the sorted quartet
        //   CF(s0,s3 | s1,s2)  = "14|23" of the sorted quartet
        //
        // We need to find which of these three sorted resolutions corresponds
        // to each of the three input CF values.

        // Build an array of the 4 indices and sort them, tracking the
        // permutation to map splits.
        struct IV {
            int idx;
            int orig;  // original position (0=t1, 1=t2, 2=t3, 3=t4)
        };
        IV arr[4] = {{i1, 0}, {i2, 1}, {i3, 2}, {i4, 3}};
        std::sort(arr, arr + 4, [](const IV& a, const IV& b) { return a.idx < b.idx; });

        // perm[original_position] = sorted_position
        int perm[4];
        for (int k = 0; k < 4; ++k) perm[arr[k].orig] = k;

        // For each of the three input CF columns, determine which sorted
        // resolution it maps to.
        //
        // Input CF12_34 is the split {orig0, orig1} | {orig2, orig3}.
        // In sorted positions, this is {perm[0], perm[1]} | {perm[2], perm[3]}.
        // The sorted resolutions are identified by the pair in the "first half":
        //   sorted 0: {s0,s1} | {s2,s3} -> pair {0,1}
        //   sorted 1: {s0,s2} | {s1,s3} -> pair {0,2}
        //   sorted 2: {s0,s3} | {s1,s2} -> pair {0,3}
        //
        // A split {a,b}|{c,d} is determined by the unordered pair {min(a,b), max(a,b)}.
        // But we need to match it to one of {0,1}, {0,2}, {0,3}.
        // Actually, a split is identified by either half. The pair that
        // contains s0 (position 0) determines the resolution:
        //   if s0 is paired with s1 -> sorted 0
        //   if s0 is paired with s2 -> sorted 1
        //   if s0 is paired with s3 -> sorted 2

        auto split_index = [&](int orig_a, int orig_b) -> int {
            // orig_a, orig_b are in the same "half" of this split.
            int sa = perm[orig_a], sb = perm[orig_b];
            // One of {sa, sb} and {the other two} forms the split.
            // Find which sorted position is paired with position 0.
            int partner_of_0 = -1;
            if (sa == 0)
                partner_of_0 = sb;
            else if (sb == 0)
                partner_of_0 = sa;
            else {
                // Position 0 is in the other half. The other half is {complement of sa,sb}.
                // complement: the two positions not in {sa, sb}.
                // Position 0 is in that complement. Its partner is the other element.
                for (int k = 0; k < 4; ++k) {
                    if (k != sa && k != sb && k != 0) {
                        partner_of_0 = k;
                        break;
                    }
                }
            }
            // partner_of_0 is 1, 2, or 3 -> maps to sorted resolution 0, 1, 2.
            return partner_of_0 - 1;
        };

        int dest_12 = split_index(0, 1);  // CF12_34: {t1,t2}|{t3,t4}
        int dest_13 = split_index(0, 2);  // CF13_24: {t1,t3}|{t2,t4}
        int dest_14 = split_index(0, 3);  // CF14_23: {t1,t4}|{t2,t3}

        std::array<int, 4> sorted_idx = {arr[0].idx, arr[1].idx, arr[2].idx, arr[3].idx};
        long long rank = rank4taxa(sorted_idx);
        size_t ri = static_cast<size_t>(rank);

        double sorted_cf[3];
        sorted_cf[dest_12] = r.cf12;
        sorted_cf[dest_13] = r.cf13;
        sorted_cf[dest_14] = r.cf14;

        cfs(ri, 0) = sorted_cf[0];
        cfs(ri, 1) = sorted_cf[1];
        cfs(ri, 2) = sorted_cf[2];
        present[ri] = true;

        if (col_ngenes >= 0) {
            ngenes_vec[ri] = r.ngenes;
            if (r.ngenes > max_ngenes) max_ngenes = r.ngenes;
        }
    }

    // Build weights: ngenes / max_ngenes for present quartets, 0 for missing.
    CsvCF result;
    result.cfs = std::move(cfs);
    result.taxa = std::move(taxa);
    if (col_ngenes >= 0 && max_ngenes > 0) {
        result.weights.resize(static_cast<size_t>(nq));
        for (size_t i = 0; i < static_cast<size_t>(nq); ++i) {
            if (present[i])
                result.weights[i] = ngenes_vec[i] / max_ngenes;
            else
                result.weights[i] = 0.0;
        }
    }
    return result;
}

bool looks_like_csv(const std::string& filename) {
    std::ifstream in(filename);
    if (!in) return false;
    std::string line;
    while (std::getline(in, line)) {
        std::string s = strip(line);
        if (s.empty()) continue;
        if (s.size() >= 3 && static_cast<unsigned char>(s[0]) == 0xEF &&
            static_cast<unsigned char>(s[1]) == 0xBB &&
            static_cast<unsigned char>(s[2]) == 0xBF) {
            s = strip(s.substr(3));
        }
        if (!s.empty() && s[0] == '(') return false;
        if (s.find(',') == std::string::npos) return false;
        std::string lower;
        lower.reserve(s.size());
        for (char c : s) lower.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
        if (lower.find("cf") != std::string::npos &&
            (lower.find("taxon") != std::string::npos || lower.find("t1") != std::string::npos ||
             lower.find("taxa") != std::string::npos)) {
            return true;
        }
        return false;
    }
    return false;
}

}  // namespace snaqpp
