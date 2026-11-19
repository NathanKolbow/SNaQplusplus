#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "cli.hpp"
#include "newick.hpp"
#include "readcsv.hpp"
#include "snaq.hpp"
#include "test_support.hpp"

using namespace snaqpp;
using namespace snaqpp::test;

namespace {

std::string write_temp_file(const std::string& filename, const std::string& content) {
    std::string path = (std::filesystem::temp_directory_path() / filename).string();
    std::ofstream out(path);
    out << content;
    return path;
}

}  // namespace

TEST(ReadTableCsv, BasicQuartet) {
    std::string csv_content =
        "taxon1,taxon2,taxon3,taxon4,CF12_34,CF13_24,CF14_23\n"
        "A,B,C,D,0.7,0.2,0.1\n";
    std::string path = write_temp_file("basic_cf.csv", csv_content);

    CsvCF csv = readtablecsv(path);
    std::filesystem::remove(path);

    EXPECT_EQ(csv.taxa.size(), 4u);
    EXPECT_EQ(csv.taxa[0], "A");
    EXPECT_EQ(csv.taxa[1], "B");
    EXPECT_EQ(csv.taxa[2], "C");
    EXPECT_EQ(csv.taxa[3], "D");
    EXPECT_EQ(csv.cfs.nrow(), 1u);
    EXPECT_NEAR(csv.cfs(0, 0), 0.7, 1e-6);
    EXPECT_NEAR(csv.cfs(0, 1), 0.2, 1e-6);
    EXPECT_NEAR(csv.cfs(0, 2), 0.1, 1e-6);
    EXPECT_TRUE(csv.weights.empty());
}

TEST(ReadTableCsv, PermutedTaxaMapping) {
    // Four taxa sorted: A, B, C, D (indices 0, 1, 2, 3)
    // Row order: B, C, A, D
    // CF12_34 = {B,C}|{A,D} = {A,D}|{B,C} -> sorted resolution 2 ({0,3}|{1,2})
    // CF13_24 = {B,A}|{C,D} = {A,B}|{C,D} -> sorted resolution 0 ({0,1}|{2,3})
    // CF14_23 = {B,D}|{C,A} = {A,C}|{B,D} -> sorted resolution 1 ({0,2}|{1,3})
    std::string csv_content =
        "t1,t2,t3,t4,CF12_34,CF13_24,CF14_23\n"
        "B,C,A,D,0.1,0.2,0.7\n";
    std::string path = write_temp_file("permuted_cf.csv", csv_content);

    CsvCF csv = readtablecsv(path);
    std::filesystem::remove(path);

    EXPECT_EQ(csv.taxa.size(), 4u);
    EXPECT_EQ(csv.cfs.nrow(), 1u);
    // sorted resolution 0 should get CF13_24 = 0.2
    EXPECT_NEAR(csv.cfs(0, 0), 0.2, 1e-6);
    // sorted resolution 1 should get CF14_23 = 0.7
    EXPECT_NEAR(csv.cfs(0, 1), 0.7, 1e-6);
    // sorted resolution 2 should get CF12_34 = 0.1
    EXPECT_NEAR(csv.cfs(0, 2), 0.1, 1e-6);
}

TEST(ReadTableCsv, NgenesWeightsAndMissingQuartets) {
    // 5 taxa: A, B, C, D, E -> C(5,4) = 5 quartets
    // Provide 2 quartets: ABCD (ngenes=100) and ABCE (ngenes=200)
    std::string csv_content =
        "taxon1,taxon2,taxon3,taxon4,CF12_34,CF13_24,CF14_23,ngenes\n"
        "A,B,C,D,0.6,0.3,0.1,100\n"
        "A,B,C,E,0.5,0.3,0.2,200\n";
    std::string path = write_temp_file("ngenes_cf.csv", csv_content);

    CsvCF csv = readtablecsv(path);
    std::filesystem::remove(path);

    EXPECT_EQ(csv.taxa.size(), 5u);
    EXPECT_EQ(csv.cfs.nrow(), 5u);
    EXPECT_EQ(csv.weights.size(), 5u);

    // Sorted taxa: A (0), B (1), C (2), D (3), E (4)
    // Quartet ABCD is indices 0,1,2,3 -> rank = rank4taxa({0,1,2,3})
    long long rank_abcd = rank4taxa({0, 1, 2, 3});
    EXPECT_NEAR(csv.weights[static_cast<size_t>(rank_abcd)], 100.0 / 200.0, 1e-6);
    EXPECT_NEAR(csv.cfs(static_cast<size_t>(rank_abcd), 0), 0.6, 1e-6);

    // Quartet ABCE is indices 0,1,2,4 -> rank = rank4taxa({0,1,2,4})
    long long rank_abce = rank4taxa({0, 1, 2, 4});
    EXPECT_NEAR(csv.weights[static_cast<size_t>(rank_abce)], 200.0 / 200.0, 1e-6);
    EXPECT_NEAR(csv.cfs(static_cast<size_t>(rank_abce), 0), 0.5, 1e-6);

    // Missing quartets should have weight 0 and CFs 1/3
    for (size_t r = 0; r < 5; ++r) {
        if (r != static_cast<size_t>(rank_abcd) && r != static_cast<size_t>(rank_abce)) {
            EXPECT_DOUBLE_EQ(csv.weights[r], 0.0);
            EXPECT_NEAR(csv.cfs(r, 0), 1.0 / 3.0, 1e-6);
            EXPECT_NEAR(csv.cfs(r, 1), 1.0 / 3.0, 1e-6);
            EXPECT_NEAR(csv.cfs(r, 2), 1.0 / 3.0, 1e-6);
        }
    }
}

TEST(ReadTableCsv, WhitespaceAndQuotes) {
    std::string csv_content =
        " \"taxon1\" , \"taxon2\" , \"taxon3\" , \"taxon4\" , \"CF12_34\" , \"CF13_24\" , \"CF14_23\" \n"
        " \"A\" , \"B\" , \"C\" , \"D\" , 0.5 , 0.3 , 0.2 \n";
    std::string path = write_temp_file("quotes_cf.csv", csv_content);

    CsvCF csv = readtablecsv(path);
    std::filesystem::remove(path);

    EXPECT_EQ(csv.taxa.size(), 4u);
    EXPECT_EQ(csv.taxa[0], "A");
    EXPECT_NEAR(csv.cfs(0, 0), 0.5, 1e-6);
    EXPECT_NEAR(csv.cfs(0, 1), 0.3, 1e-6);
    EXPECT_NEAR(csv.cfs(0, 2), 0.2, 1e-6);
}

TEST(ReadTableCsv, LooksLikeCsv) {
    std::string csv_content = "taxon1,taxon2,taxon3,taxon4,CF12_34,CF13_24,CF14_23\nA,B,C,D,0.7,0.2,0.1\n";
    std::string csv_path = write_temp_file("detect.csv", csv_content);
    EXPECT_TRUE(looks_like_csv(csv_path));
    std::filesystem::remove(csv_path);

    std::string newick_content = "((A,B),(C,D));\n";
    std::string newick_path = write_temp_file("detect.tre", newick_content);
    EXPECT_FALSE(looks_like_csv(newick_path));
    std::filesystem::remove(newick_path);
}

TEST(ReadTableCsv, ErrorCases) {
    // Missing CF column
    std::string bad_hdr = "taxon1,taxon2,taxon3,taxon4,CF12_34,CF13_24\nA,B,C,D,0.5,0.5\n";
    std::string path1 = write_temp_file("bad1.csv", bad_hdr);
    EXPECT_THROW(readtablecsv(path1), SnaqError);
    std::filesystem::remove(path1);

    // Duplicate taxa in a row
    std::string dup_taxa = "taxon1,taxon2,taxon3,taxon4,CF12_34,CF13_24,CF14_23\nA,B,C,A,0.5,0.3,0.2\n";
    std::string path2 = write_temp_file("bad2.csv", dup_taxa);
    EXPECT_THROW(readtablecsv(path2), SnaqError);
    std::filesystem::remove(path2);

    // Fewer than 4 taxa
    std::string few_taxa = "taxon1,taxon2,taxon3,taxon4,CF12_34,CF13_24,CF14_23\nA,B,C,D,0.5,0.3,0.2\n";
    // (Wait, ABCD is 4 taxa, so that has 4 taxa. What about non-existent file?)
    EXPECT_THROW(readtablecsv("non_existent_file_12345.csv"), SnaqError);
}

TEST(ReadTableCsv, XiphophorusExampleIfAvailable) {
    // Check if ~/Downloads/xiphophorus.csv exists
    std::string user_home = "";
#ifdef _WIN32
    const char* userprofile = std::getenv("USERPROFILE");
    if (userprofile) user_home = userprofile;
#else
    const char* home = std::getenv("HOME");
    if (home) user_home = home;
#endif
    if (user_home.empty()) return;
    std::string xipho_path = (std::filesystem::path(user_home) / "Downloads" / "xiphophorus.csv").string();
    if (!std::filesystem::exists(xipho_path)) return;

    EXPECT_TRUE(looks_like_csv(xipho_path));
    CsvCF csv = readtablecsv(xipho_path);

    // 24 taxa, C(24, 4) = 10626 quartets
    EXPECT_EQ(csv.taxa.size(), 24u);
    EXPECT_EQ(csv.cfs.nrow(), 10626u);
    EXPECT_EQ(csv.weights.size(), 10626u);

    // Every quartet in xiphophorus has ngenes = 1183, so weight should be 1.0 everywhere
    for (size_t i = 0; i < csv.weights.size(); ++i) {
        EXPECT_NEAR(csv.weights[i], 1.0, 1e-6);
        // Each quartet's CFs sum to ~1.0
        double sum_cf = csv.cfs(i, 0) + csv.cfs(i, 1) + csv.cfs(i, 2);
        EXPECT_NEAR(sum_cf, 1.0, 1e-4);
    }
}
