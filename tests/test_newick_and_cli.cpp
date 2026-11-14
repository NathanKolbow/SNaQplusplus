// Tests specific to the C++ port: Newick I/O, Julia number formatting,
// quartet ranking and command-line parsing.
#include "cli.hpp"
#include "quartets.hpp"
#include "test_support.hpp"

using namespace csnaq;
using namespace csnaq::test;

TEST(Newick, RoundTrip) {
    for (const char* s : {"((a:1.0,b:2.5):0.5,(c:1.0,d:1.0):0.25,e:3.0);",
                          "(((a:1.0,(b:1.0)#H1:0.5::0.7):1.0,(c:1.0,#H1:0.2::0.3):1.0):1.0,d:1.0);",
                          "((a,b)i1,((c,#H1)i6,(d,((e,f)i5)#H1)i3)i2);"}) {
        Network net = nw(s);
        EXPECT_EQ(writenewick(net), s);
    }
}

TEST(Newick, RootEdgeAndComments) {
    Network net = nw("((a:1,b:1)[&comment]:1,c:2):5;");
    EXPECT_EQ(net.numtaxa, 3);
    EXPECT_EQ(writenewick(net), "((a:1.0,b:1.0):1.0,c:2.0);");
}

TEST(Newick, Errors) {
    EXPECT_THROW(nw("((a,b),c)"), SnaqError);
    EXPECT_THROW(nw("((a,b),c));"), SnaqError);
}

TEST(Newick, JuliaFloatFormat) {
    EXPECT_EQ(julia_float(0.0), "0.0");
    EXPECT_EQ(julia_float(1e-4), "0.0001");
    EXPECT_EQ(julia_float(1e-5), "1.0e-5");
    EXPECT_EQ(julia_float(123456.0), "123456.0");
    EXPECT_EQ(julia_float(1234567.0), "1.234567e6");
    EXPECT_EQ(julia_float(0.25995127997690437), "0.25995127997690437");
    EXPECT_EQ(julia_float(-0.2835919482810475), "-0.2835919482810475");
}

TEST(Quartets, RankUnrank) {
    for (int n : {4, 5, 9, 30}) {
        long long nq = nchoose4(n);
        std::array<int, 4> prev{-1, -1, -1, -1};
        for (long long r = 0; r < nq; ++r) {
            auto q = unrank4taxa(n, r);
            EXPECT_TRUE(q[0] < q[1] && q[1] < q[2] && q[2] < q[3] && q[3] < n);
            EXPECT_EQ(rank4taxa(q), r);
            // colexicographic order
            if (r > 0) EXPECT_TRUE(std::lexicographical_compare(prev.rbegin(), prev.rend(), q.rbegin(), q.rend()));
            prev = q;
        }
    }
}

TEST(CLI, ParsesAllSnaqKeywords) {
    std::vector<std::string> args{"-i", "g.tre", "--start-tree", "(a,b,(c,d));", "-o", "out/x", "--hmax", "3",
                                  "--runs", "7", "--Nfail", "11", "--outgroup", "a", "--seed", "-5", "--probST",
                                  "0.5", "--updateBL", "false", "--probQR", "0.25", "--propQuartets", "0.5",
                                  "--propQuartetsFinal=0.75", "--qinfTest", "--qtolAbs", "0.001", "--ftolRel",
                                  "1e-6", "--ftolAbs", "1e-5", "--xtolRel", "1e-2", "--xtolAbs", "1e-3", "--rho",
                                  "0.3", "--verbose", "--max_level", "2", "--galled_network", "--maxeval", "100",
                                  "--liktolAbs", "1e-7", "--liktolRel", "1e-3", "--opt_maxeval", "40",
                                  "--logprefix", "lp", "--lazy", "--threads", "2"};
    CliArgs a = parse_cli(args);
    const SnaqOptions& o = a.opts;
    EXPECT_EQ(a.input, "g.tre");
    EXPECT_EQ(a.start_tree, "(a,b,(c,d));");
    EXPECT_EQ(o.filename, "out/x");
    EXPECT_EQ(o.hmax, 3);
    EXPECT_EQ(o.runs, 7);
    EXPECT_EQ(o.Nfail, 11);
    EXPECT_EQ(o.outgroup, "a");
    EXPECT_EQ(o.seed, -5);
    EXPECT_TRUE(o.seed_given);
    EXPECT_DOUBLE_EQ(o.probST, 0.5);
    EXPECT_FALSE(o.updateBL);
    EXPECT_DOUBLE_EQ(o.probQR, 0.25);
    EXPECT_DOUBLE_EQ(o.propQuartets, 0.5);
    EXPECT_DOUBLE_EQ(o.propQuartetsFinal, 0.75);
    EXPECT_TRUE(o.qinfTest);
    EXPECT_DOUBLE_EQ(o.qtolAbs, 0.001);
    EXPECT_DOUBLE_EQ(o.ftolRel, 1e-6);
    EXPECT_DOUBLE_EQ(o.ftolAbs, 1e-5);
    EXPECT_DOUBLE_EQ(o.xtolRel, 1e-2);
    EXPECT_DOUBLE_EQ(o.xtolAbs, 1e-3);
    EXPECT_DOUBLE_EQ(o.rho, 0.3);
    EXPECT_TRUE(o.verbose);
    EXPECT_EQ(o.restrictions.max_level, 2);
    EXPECT_TRUE(o.restrictions.galled_network);
    EXPECT_EQ(o.maxeval, 100);
    EXPECT_DOUBLE_EQ(o.liktolAbs, 1e-7);
    EXPECT_DOUBLE_EQ(o.liktolRel, 1e-3);
    EXPECT_EQ(o.opt_maxeval, 40);
    EXPECT_EQ(o.logprefix, "lp");
    EXPECT_TRUE(o.lazy);
    EXPECT_EQ(o.threads, 2);
}

TEST(CLI, LazyIsAutomaticUnlessGiven) {
    // not given: lazy is decided from propQuartets and propQuartetsFinal (snaq)
    CliArgs a = parse_cli({"-i", "g.tre", "--start-tree", "t.tre"});
    EXPECT_FALSE(a.opts.lazy);
    EXPECT_TRUE(a.opts.autoLazy);
    CliArgs b = parse_cli({"-i", "g.tre", "--start-tree", "t.tre", "--lazy"});
    EXPECT_TRUE(b.opts.lazy);
    CliArgs c = parse_cli({"-i", "g.tre", "--start-tree", "t.tre", "--lazy", "false"});
    EXPECT_FALSE(c.opts.lazy);
    EXPECT_FALSE(c.opts.autoLazy);
    CliArgs d = parse_cli({"-i", "g.tre", "--start-tree", "t.tre", "--lazy=true", "--propQuartets", "0.5"});
    EXPECT_TRUE(d.opts.lazy);
}

TEST(CLI, DefaultsFollowSnaqCode) {
    CliArgs a = parse_cli({"-i", "g.tre", "--start-tree", "t.tre"});
    const SnaqOptions& o = a.opts;
    EXPECT_EQ(o.hmax, 1);
    EXPECT_EQ(o.runs, 100);
    EXPECT_EQ(o.Nfail, 50);
    EXPECT_EQ(o.filename, "snaq");
    EXPECT_EQ(o.outgroup, "none");
    EXPECT_DOUBLE_EQ(o.ftolRel, 1e-8);
    EXPECT_DOUBLE_EQ(o.xtolAbs, 1e-8);
    EXPECT_DOUBLE_EQ(o.probST, 0.9);
    EXPECT_TRUE(o.updateBL);
    EXPECT_FALSE(o.seed_given);
    EXPECT_FALSE(o.sweep);
    EXPECT_TRUE(parse_cli({"-i", "g.tre", "--start-tree", "t.tre", "--sweep"}).opts.sweep);
    EXPECT_FALSE(parse_cli({"-i", "g.tre", "--start-tree", "t.tre", "--sweep", "false"}).opts.sweep);
}

TEST(CLI, UsageErrors) {
    EXPECT_THROW(parse_cli({"--start-tree", "t.tre"}), UsageError);
    EXPECT_THROW(parse_cli({"-i", "g.tre"}), UsageError);
    EXPECT_THROW(parse_cli({"-i", "g.tre", "--start-tree", "t", "--bogus"}), UsageError);
    EXPECT_THROW(parse_cli({"-i", "g.tre", "--start-tree", "t", "--hmax", "x"}), UsageError);
}
