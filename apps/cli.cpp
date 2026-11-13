#include "cli.hpp"

#include <functional>
#include <map>
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
        throw ArgumentError("invalid value for " + flag + ": '" + v + "' (expected a number)");
    }
}

long long parse_int(const std::string& flag, const std::string& v) {
    try {
        size_t idx = 0;
        long long d = std::stoll(v, &idx);
        if (idx != v.size()) throw std::invalid_argument(v);
        return d;
    } catch (...) {
        throw ArgumentError("invalid value for " + flag + ": '" + v + "' (expected an integer)");
    }
}

bool is_bool_word(const std::string& v) { return v == "true" || v == "false"; }

}  // namespace

std::string usage() {
    return R"(snaq: phylogenetic network inference by maximum pseudolikelihood (C++ port of SNaQ.jl's snaq!)

usage:
  snaq -i <gene trees> -o <output root name> --start-tree <newick or file> [options]

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

advanced (keyword arguments snaq! passes to each search):
  --maxeval INT             maximum number of search iterations (default: 100000000)
  --liktolAbs FLOAT         absolute score improvement to accept a move (default: 1e-8)
  --liktolRel FLOAT         relative score improvement to accept a move (default: 1e-4)
  --opt_maxeval INT         max. evaluations when fitting a proposed network (default: max(30, #taxa))
  --logprefix PATH          write a detailed per-run log to PATH<run number>

other:
  --lazy                    compute observed CFs on demand (like DataCF(trees; lazy=true));
                            requires --propQuartets < 1
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
    integer("--maxeval", &o.maxeval);
    num("--liktolAbs", &o.liktolAbs);
    num("--liktolRel", &o.liktolRel);
    integer("--opt_maxeval", &o.opt_maxeval);
    str("--logprefix", &o.logprefix);
    boolean("--lazy", &o.lazy);
    integer("--threads", &o.threads);
    // hyphenated aliases for the restriction / advanced flags
    for (const char* name : {"max_level", "galled_tree", "galled_network", "rooted_tree_child", "weakly_tree_child",
                             "strongly_tree_child", "opt_maxeval"}) {
        std::string alias = std::string("--") + name;
        for (auto& c : alias)
            if (c == '_') c = '-';
        specs[alias] = specs[std::string("--") + name];
    }

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
        if (it == specs.end()) throw ArgumentError("unknown option: " + args[k] + " (see snaq --help)");
        const Spec& sp = it->second;
        if (sp.is_bool) {
            if (!has_value) {
                if (k + 1 < args.size() && is_bool_word(args[k + 1])) {
                    value = args[++k];
                } else {
                    value = "true";
                }
            } else if (!is_bool_word(value)) {
                throw ArgumentError("invalid value for " + tok + ": '" + value + "' (expected true or false)");
            }
        } else if (!has_value) {
            if (k + 1 >= args.size()) throw ArgumentError("missing value for " + tok);
            value = args[++k];
        }
        sp.set(tok, value);
    }
    if (a.help || a.version) return a;
    if (a.input.empty()) throw ArgumentError("missing required option -i/--input (gene trees)");
    if (a.start_tree.empty()) throw ArgumentError("missing required option --start-tree");
    if (o.threads < 1) throw ArgumentError("--threads must be ≥ 1");
    return a;
}

}  // namespace csnaq
