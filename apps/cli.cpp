#include "cli.hpp"

#include <functional>
#include <map>
#include <set>
#include <thread>

namespace csnaq {

namespace {

double parse_double(const std::string& flag, const std::string& v) {
    try {
        size_t idx = 0;
        double d = std::stod(v, &idx);
        if (idx != v.size()) throw std::invalid_argument(v);
        return d;
    } catch (...) {
        throw UsageError("invalid value for " + flag + ": '" + v + "' (expected a number)");
    }
}

long long parse_int(const std::string& flag, const std::string& v) {
    try {
        size_t idx = 0;
        long long d = std::stoll(v, &idx);
        if (idx != v.size()) throw std::invalid_argument(v);
        return d;
    } catch (...) {
        throw UsageError("invalid value for " + flag + ": '" + v + "' (expected an integer)");
    }
}

bool is_bool_word(const std::string& v) { return v == "true" || v == "false"; }

}  // namespace

std::string short_usage() {
    return R"(USAGE: snaq -i <gene trees file> --start-tree <newick string or file> [-o <output root name>] [options]

Infers a phylogenetic network from gene trees (C++ port of SNaQ.jl's snaq!).

Common options:
  -i, --input FILE          gene trees, one Newick tree per line (required)
  --start-tree TREE         starting tree: Newick string or file (required)
  -o, --filename PATH       root name of the output files (default: snaq)
  --hmax INT                maximum number of hybridizations (default: 1)
  --runs INT                number of independent runs (default: 100)
  --Nfail INT               failed proposals in a row before a run stops (default: 50)
  --outgroup NAME           root the networks at this taxon
  --seed INT                random seed
  --threads INT             number of threads (default: all cores)

Example:
  snaq -i genetrees.tre --start-tree astral.tre -o net1 --hmax 1 --runs 10

Optimize only the branch lengths and γ of a fixed network:
  snaq --optimize net.tre -i genetrees.tre -o net_optimized.tre [--propQuartets FLOAT]

Run `snaq --help` for all options.
)";
}

std::string usage() {
    return R"(snaq: phylogenetic network inference by maximum pseudolikelihood (C++ port of SNaQ.jl's snaq!)

usage:
  snaq -i <gene trees> -o <output root name> --start-tree <newick or file> [options]
  snaq --optimize <network> -i <gene trees> [-o <output file>] [optimization options]

required:
  -i, --input FILE          gene trees, one Newick tree per line
  --start-tree TREE         starting topology: a Newick string (e.g. "((a,b),c,d);") or a file
                            with one or more Newick lines (one per run, or one for all runs)

output:
  -o, --filename PATH       root name for output files .out, .networks, .log and _runs/
                            (default: snaq). An empty string ("") writes no files.

snaq! options (same names and defaults as SNaQ.jl's snaq!):
  --hmax INT                maximum number of hybridizations (default: 1)
  --runs INT                number of independent runs (default: 100)
  --Nfail INT               consecutive rejected proposals before a run stops (default: 50)
  --outgroup NAME           outgroup taxon to root the networks at (default: none)
  --seed INT                random seed (default: random)
  --probST FLOAT            probability of perturbing the start topology by an NNI (default: 0.9)
  --updateBL [BOOL]         optimize parameters of the start topology first (default: true)
  --probQR FLOAT            probability of quartet-guided move proposals (default: 0.0)
  --propQuartets FLOAT      proportion of quartets used during the search (default: 1.0)
  --propQuartetsFinal FLOAT proportion of quartets for the final re-optimization (default: 1.0)
  --qinfTest [BOOL]         discard uninformative quartets (default: false)
  --qtolAbs FLOAT           tolerance of the uninformative-quartet test (default: 1e-4)
  --ftolRel FLOAT           relative score tolerance of parameter optimization (default: 1e-8)
  --ftolAbs FLOAT           absolute score tolerance (default: 1e-8)
  --xtolRel FLOAT           relative parameter tolerance (default: 1e-8)
  --xtolAbs FLOAT           absolute parameter tolerance (default: 1e-8)
  --rho FLOAT               inheritance correlation ρ in [0, 1] (default: 0.0)
  --verbose [BOOL]          print progress of each run (default: false)

restrictions (as in SNaQ.jl's restrictionset; default: none):
  --max_level INT           maximum network level
  --galled_tree             galled trees (level-1 networks)
  --galled_network          galled networks
  --rooted_tree_child       rooted tree-child networks
  --weakly_tree_child       weakly tree-child networks
  --strongly_tree_child     strongly tree-child networks
  --tcg_identifiable        networks passing SNaQ.jl's tcgidentifiable (identifiable
                            topology and branch parameters; trees always pass)

advanced (keyword arguments snaq! passes to each search):
  --maxeval INT             maximum number of search iterations (default: 100000000)
  --liktolAbs FLOAT         absolute score improvement to accept a move (default: 1e-8)
  --liktolRel FLOAT         relative score improvement to accept a move (default: 1e-4)
  --opt_maxeval INT         max. evaluations when fitting a proposed network
                            (default: 200 with newton, max(30, #taxa) with lbfgs)
  --logprefix PATH          write a detailed per-run log to PATH<run number>

pseudolikelihood (not in SNaQ.jl):
  --geneWeights [BOOL]      weigh each quartet by the proportion of gene trees its observed CFs
                            come from, instead of weighing all quartets equally (default: false)

parameter optimization (not in SNaQ.jl):
  --optimizer NAME          "newton": projected Newton method, run to convergence (default);
                            "lbfgs": L-BFGS limited to --opt_maxeval evaluations, as SNaQ.jl

hybrid additions (not in SNaQ.jl):
  --guideHybrids [BOOL]     draw half of the proposed hybrid additions between edges where the
                            data suggest a reticulation (default: true)
  --minAddHybrid INT        with fewer than hmax hybridizations, a run does not stop before this
                            many hybrid additions were evaluated since its last accepted move
                            (default: Nfail; 0: no such requirement)

local search (not in SNaQ.jl):
  --sweep [BOOL]            before a run stops, evaluate every local move of its network (the
                            NNIs and the local moves of its reticulations); the run goes on if
                            one is accepted (default: false)

search shortcuts (not in SNaQ.jl; each can be turned off with "false"):
  --boundPrune [BOOL]       reject a proposed network without fitting it when an upper bound of
                            its score shows it cannot be accepted (default: true)
  --dedup [BOOL]            do not fit again a topology that was already fitted (default: true)
  --dedupCounts [BOOL]      count such repeated topologies as failed proposals, towards --Nfail
                            (default: false: another proposal is drawn instead)

parameter optimization of a fixed network (not in SNaQ.jl):
  --optimize NET            do not search: optimize only the branch lengths and inheritance
                            probabilities γ of the topology NET (a Newick string or a file with
                            one Newick line). -o is then the file the optimized network is
                            written to (default: none, the network is only printed).
                            Other options allowed with --optimize: --propQuartets (proportion
                            of quartets used, sampled with --seed), --optimizer, --opt_maxeval,
                            --ftolRel, --ftolAbs, --xtolRel, --xtolAbs, --rho, --geneWeights,
                            --outgroup, --threads

other:
  --lazy [BOOL]             compute observed CFs on demand (like DataCF(trees; lazy=true));
                            requires --propQuartets < 1. Default: on when --propQuartets and
                            --propQuartetsFinal are both below 1, off otherwise
  --threads INT             number of threads (default: all available cores)
  -h, --help                show this help
  --version                 show the version
)";
}

CliArgs parse_cli(const std::vector<std::string>& args) {
    CliArgs a;
    SnaqOptions& o = a.opts;
    unsigned hw = std::thread::hardware_concurrency();
    o.threads = hw > 0 ? static_cast<int>(hw) : 1;

    using Setter = std::function<void(const std::string& flag, const std::string& value)>;
    struct Spec {
        bool is_bool;
        Setter set;
    };
    std::map<std::string, Spec> specs;
    auto num = [&](const std::string& name, double* dst) {
        specs[name] = {false, [dst](const std::string& f, const std::string& v) { *dst = parse_double(f, v); }};
    };
    auto integer = [&](const std::string& name, auto* dst) {
        specs[name] = {false, [dst](const std::string& f, const std::string& v) {
                           *dst = static_cast<std::remove_pointer_t<decltype(dst)>>(parse_int(f, v));
                       }};
    };
    auto str = [&](const std::string& name, std::string* dst) {
        specs[name] = {false, [dst](const std::string&, const std::string& v) { *dst = v; }};
    };
    auto boolean = [&](const std::string& name, bool* dst) {
        specs[name] = {true, [dst](const std::string&, const std::string& v) { *dst = (v == "true"); }};
    };

    str("-i", &a.input);
    str("--input", &a.input);
    str("--start-tree", &a.start_tree);
    str("--optimize", &a.optimize);
    str("-o", &o.filename);
    str("--filename", &o.filename);
    integer("--hmax", &o.hmax);
    integer("--runs", &o.runs);
    integer("--Nfail", &o.Nfail);
    str("--outgroup", &o.outgroup);
    specs["--seed"] = {false, [&o](const std::string& f, const std::string& v) {
                           o.seed = parse_int(f, v);
                           o.seed_given = true;
                       }};
    num("--probST", &o.probST);
    boolean("--updateBL", &o.updateBL);
    num("--probQR", &o.probQR);
    num("--propQuartets", &o.propQuartets);
    num("--propQuartetsFinal", &o.propQuartetsFinal);
    boolean("--qinfTest", &o.qinfTest);
    num("--qtolAbs", &o.qtolAbs);
    num("--ftolRel", &o.ftolRel);
    num("--ftolAbs", &o.ftolAbs);
    num("--xtolRel", &o.xtolRel);
    num("--xtolAbs", &o.xtolAbs);
    num("--rho", &o.rho);
    boolean("--verbose", &o.verbose);
    integer("--max_level", &o.restrictions.max_level);
    boolean("--galled_tree", &o.restrictions.galled_tree);
    boolean("--galled_network", &o.restrictions.galled_network);
    boolean("--rooted_tree_child", &o.restrictions.rooted_tree_child);
    boolean("--weakly_tree_child", &o.restrictions.weakly_tree_child);
    boolean("--strongly_tree_child", &o.restrictions.strongly_tree_child);
    boolean("--tcg_identifiable", &o.restrictions.tcg_identifiable);
    integer("--maxeval", &o.maxeval);
    num("--liktolAbs", &o.liktolAbs);
    num("--liktolRel", &o.liktolRel);
    integer("--opt_maxeval", &o.opt_maxeval);
    str("--logprefix", &o.logprefix);
    boolean("--geneWeights", &o.geneWeights);
    str("--optimizer", &o.optimizer);
    boolean("--guideHybrids", &o.guideHybrids);
    integer("--minAddHybrid", &o.minAddHybrid);
    boolean("--sweep", &o.sweep);
    boolean("--boundPrune", &o.boundPrune);
    boolean("--dedup", &o.dedup);
    boolean("--dedupCounts", &o.dedupCounts);
    // --lazy false also turns off the automatic use of on-demand CFs
    specs["--lazy"] = {true, [&o](const std::string&, const std::string& v) {
                           o.lazy = (v == "true");
                           o.autoLazy = false;
                       }};
    integer("--threads", &o.threads);
    // hyphenated aliases for the restriction / advanced flags
    for (const char* name : {"max_level", "galled_tree", "galled_network", "rooted_tree_child", "weakly_tree_child",
                             "strongly_tree_child", "tcg_identifiable", "opt_maxeval"}) {
        std::string alias = std::string("--") + name;
        for (auto& c : alias)
            if (c == '_') c = '-';
        specs[alias] = specs[std::string("--") + name];
    }

    std::set<std::string> seen;
    for (size_t k = 0; k < args.size(); ++k) {
        std::string tok = args[k];
        if (tok == "-h" || tok == "--help") {
            a.help = true;
            continue;
        }
        if (tok == "--version") {
            a.version = true;
            continue;
        }
        std::string value;
        bool has_value = false;
        auto eq = tok.find('=');
        if (tok.rfind("--", 0) == 0 && eq != std::string::npos) {
            value = tok.substr(eq + 1);
            tok = tok.substr(0, eq);
            has_value = true;
        }
        auto it = specs.find(tok);
        if (it == specs.end()) throw UsageError("unknown option: " + args[k] + " (see snaq --help)");
        const Spec& sp = it->second;
        if (sp.is_bool) {
            if (!has_value) {
                if (k + 1 < args.size() && is_bool_word(args[k + 1])) {
                    value = args[++k];
                } else {
                    value = "true";
                }
            } else if (!is_bool_word(value)) {
                throw UsageError("invalid value for " + tok + ": '" + value + "' (expected true or false)");
            }
        } else if (!has_value) {
            if (k + 1 >= args.size()) throw UsageError("missing value for " + tok);
            value = args[++k];
        }
        sp.set(tok, value);
        seen.insert(tok);
    }
    if (a.help || a.version) return a;
    if (!a.optimize.empty()) {
        static const std::set<std::string> allowed = {
            "-i", "--input", "-o", "--filename", "--optimize", "--propQuartets", "--seed", "--optimizer",
            "--opt_maxeval", "--opt-maxeval", "--ftolRel", "--ftolAbs", "--xtolRel", "--xtolAbs", "--rho",
            "--geneWeights", "--outgroup", "--threads"};
        for (const auto& f : seen)
            if (!allowed.count(f)) throw UsageError("option " + f + " does not apply to --optimize");
        if (a.input.empty()) throw UsageError("missing required option -i/--input (gene trees)");
        if (o.threads < 1) throw UsageError("--threads must be ≥ 1");
        a.output = (seen.count("-o") || seen.count("--filename")) ? o.filename : "";
        return a;
    }
    if (a.input.empty()) throw UsageError("missing required option -i/--input (gene trees)");
    if (a.start_tree.empty()) throw UsageError("missing required option --start-tree");
    if (o.threads < 1) throw UsageError("--threads must be ≥ 1");
    return a;
}

}  // namespace csnaq
