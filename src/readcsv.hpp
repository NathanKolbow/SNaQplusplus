// Reading quartet concordance factors from a CSV table.
#pragma once

#include <string>
#include <vector>

#include "likelihood.hpp"

namespace snaqpp {

// A CF table read from a CSV file: the CFs of all C(n,4) quartets of `taxa`,
// stored in SNaQ's quartet order (colexicographic rank of sorted 0-based taxon
// indices). `weights[rank]` is the gene count of the quartet divided by the
// maximum gene count across all quartets, or empty if the CSV has no `ngenes`
// column. Missing quartets (not in the CSV) get CFs of (1/3, 1/3, 1/3) and
// weight 0.
struct CsvCF {
    CFMatrix cfs;
    std::vector<std::string> taxa;
    std::vector<double> weights;
};

// Reads a CSV file with columns:
//   taxon1, taxon2, taxon3, taxon4, CF12_34, CF13_24, CF14_23 [, ngenes]
// The header is required and identifies the columns (the order of the columns
// does not matter). `ngenes` is optional: when present, quartets are weighted
// by ngenes / max(ngenes).
//
// Each row gives the observed CFs of one quartet. The three CF columns
// correspond to the three resolutions of the four taxa:
//   CF12_34  = CF of the split {taxon1, taxon2} | {taxon3, taxon4}
//   CF13_24  = CF of the split {taxon1, taxon3} | {taxon2, taxon4}
//   CF14_23  = CF of the split {taxon1, taxon4} | {taxon2, taxon3}
//
// The quartets do not need to cover every 4-subset of the taxa that appear,
// and the same quartet may appear more than once (the last occurrence wins).
CsvCF readtablecsv(const std::string& filename);

// Returns true if the file looks like a CF CSV table: the first non-empty line
// contains a comma and the substring "CF".
bool looks_like_csv(const std::string& filename);

}  // namespace snaqpp
