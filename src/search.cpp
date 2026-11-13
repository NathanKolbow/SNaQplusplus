#include "search.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <numeric>
#include <set>
#include <thread>

#include "classes.hpp"
#include "equations.hpp"
#include "moves.hpp"
#include "network_manip.hpp"
#include "newick.hpp"
#include "optimize.hpp"
#include "signature.hpp"
#include "snaq_network.hpp"

namespace csnaq {

namespace fs = std::filesystem;

static std::mutex g_stdout_mutex;

// ------------------------------------------------------------------------
// CFSource
// ------------------------------------------------------------------------

CFSource CFSource::full(CFMatrix m, std::vector<std::string> taxa, std::vector<double> weights) {
    CFSource s;
    if (!weights.empty() && weights.size() != m.nrow()) throw SnaqError("one weight is needed per quartet");
    s.fullweights_ = std::move(weights);
    s.full_ = std::move(m);
    s.taxa_ = std::move(taxa);
    return s;
}

CFSource CFSource::lazy(const ObservedCFs* obs) {
    CFSource s;
    s.lazy_ = obs;
    return s;
}

long long CFSource::nrows() const { return lazy_ ? lazy_->nquartets() : static_cast<long long>(full_.nrow()); }

const std::vector<std::string>& CFSource::taxa() const { return lazy_ ? lazy_->taxa() : taxa_; }

CFMatrix CFSource::rows(const std::vector<long long>& ranks, int nthreads) const {
    CFMatrix out;
    out.v.resize(3 * ranks.size());
    if (lazy_) {
        lazy_->compute_rows(ranks, out.v.data(), nthreads);
    } else {
        for (size_t j = 0; j < ranks.size(); ++j)
            for (int k = 0; k < 3; ++k) out(j, k) = full_(static_cast<size_t>(ranks[j]), k);
    }
    return out;
}

std::vector<double> CFSource::weights(const std::vector<long long>& ranks) const {
    std::vector<double> w;
    if (lazy_) {
        w.resize(ranks.size());
        for (size_t j = 0; j < ranks.size(); ++j) w[j] = lazy_->genefraction(ranks[j]);
    } else if (!fullweights_.empty()) {
        w.resize(ranks.size());
        for (size_t j = 0; j < ranks.size(); ++j) w[j] = fullweights_[static_cast<size_t>(ranks[j])];
    }
    return w;
}

const CFMatrix& CFSource::matrix() const {
    if (lazy_) throw SnaqError("the observed CFs of every quartet are not available from a lazy DataCF");
    return full_;
}

// ------------------------------------------------------------------------
// logging (search_logging.jl)
// ------------------------------------------------------------------------

void logmessage(const std::string& filename, const std::string& msg) {
    if (filename.empty()) return;
    std::ofstream f(filename + ".log", std::ios::app);
    f << msg << "\n";
}

static void logtext(const std::string& logfile, const std::string& msg) {
    if (logfile.empty()) return;
    std::ofstream f(logfile, std::ios::app);
    f << "[" << current_timestamp() << "] " << msg << "\n";
}

static size_t utf8_length(const std::string& s) {
    size_t n = 0;
    for (unsigned char c : s)
        if ((c & 0xC0) != 0x80) ++n;
    return n;
}

static std::string julia_number(double x) {
    // Julia prints whole-valued Float64 with ".0" and Int without; values here
    // are Float64 except `maximum(...; init=0)` of an empty vector (Int 0).
    return julia_float(x);
}

struct MoveStats {
    std::map<std::string, long> proposed, accepted;
    std::map<std::string, std::vector<double>> dPL;
};

static void logmoves(const std::string& logfile, const MoveStats& st) {
    if (logfile.empty() || st.proposed.empty()) return;
    size_t min_width = 0;
    for (const auto& kv : st.proposed) {
        const std::string& k = kv.first;
        min_width = std::max({min_width, utf8_length(k), std::to_string(kv.second).size(),
                              std::to_string(st.accepted.at(k)).size()});
    }
    min_width += 2;
    auto expand = [&](const std::string& s) {
        std::string r = s;
        for (size_t j = utf8_length(s); j < min_width; ++j) r += " ";
        return r;
    };
    std::string msg = std::string(12, '-') + "MOVE ACCEPTANCE RATES" + std::string(12, '-') + "\n";
    msg += expand("\tmove:");
    for (const auto& kv : st.proposed) msg += "| " + expand(kv.first);
    msg += " |\n" + expand("\tproposals:");
    for (const auto& kv : st.proposed) msg += "| " + expand(std::to_string(kv.second));
    msg += " |\n" + expand("\taccepted:");
    for (const auto& kv : st.proposed) msg += "| " + expand(std::to_string(st.accepted.at(kv.first)));
    msg += " |\n" + expand("\taccept %:");
    for (const auto& kv : st.proposed)
        msg += "| " + expand(julia_number(julia_round(100.0 * static_cast<double>(st.accepted.at(kv.first)) /
                                                          static_cast<double>(kv.second),
                                                      0)));
    msg += " |\n" + expand("\tmean SUCC ΔPL:");
    for (const auto& kv : st.proposed) {
        const auto& v = st.dPL.at(kv.first);
        double m = v.empty() ? std::numeric_limits<double>::quiet_NaN()
                             : std::accumulate(v.begin(), v.end(), 0.0) / static_cast<double>(v.size());
        msg += "| " + expand(julia_number(julia_round(m, 2)));
    }
    msg += " |\n" + expand("\tmax SUCC ΔPL:");
    for (const auto& kv : st.proposed) {
        const auto& v = st.dPL.at(kv.first);
        if (v.empty()) {
            msg += "| " + expand("0");  // maximum(...; init=0) of an empty vector is the Int 0
        } else {
            double m = 0.0;
            for (double x : v) m = std::max(m, x);
            msg += "| " + expand(julia_number(julia_round(m, 2)));
        }
    }
    msg += " |\n" + expand("\tmean ALL ΔPL:");
    for (const auto& kv : st.proposed) {
        const auto& v = st.dPL.at(kv.first);
        double s = std::accumulate(v.begin(), v.end(), 0.0);
        msg += "| " + expand(julia_number(julia_round(s / static_cast<double>(kv.second), 2)));
    }
    msg += " |\n\n";
    logtext(logfile, msg);
}

// ------------------------------------------------------------------------
// search
// ------------------------------------------------------------------------

static FitMethod fit_method(const SnaqOptions& o) {
    if (o.optimizer == "newton") return FitMethod::newton;
    if (o.optimizer == "lbfgs") return FitMethod::lbfgs;
    throw SnaqError("optimizer must be \"newton\" or \"lbfgs\" (optimizer = " + o.optimizer + ")");
}

static FitOptions search_fit_options(const SnaqOptions& o, int maxeval) {
    FitOptions f;
    f.method = fit_method(o);
    f.maxeval = maxeval;
    f.ftolRel = o.ftolRel;
    f.ftolAbs = o.ftolAbs;
    f.xtolRel = o.xtolRel;
    f.xtolAbs = o.xtolAbs;
    return f;
}

static std::vector<long long> allranks(long long n) {
    std::vector<long long> r(static_cast<size_t>(n));
    std::iota(r.begin(), r.end(), 0LL);
    return r;
}

static std::string round5(double x) { return julia_float(julia_round(x, 5)); }
static std::string round3(double x) { return julia_float(julia_round(x, 3)); }
static std::string round6(double x) { return julia_float(julia_round(x, 6)); }

Network search(Network N, const CFSource& q, int hmax, const SnaqOptions& o, long long seed,
               double propQuartetsFinal, const std::string& logfile, const std::string& filename, int nthreads,
               SearchShared* shared) {
    const double rho = o.rho;
    const double propQuartets = o.propQuartets;
    const long long maxequivPLs = o.Nfail;
    // (with the Newton method, -1 lets fits run to convergence: see FitOptions)
    const int opt_maxeval =
        o.opt_maxeval >= 0 ? o.opt_maxeval : (fit_method(o) == FitMethod::newton ? -1 : std::max(30, N.numtaxa));
    if (!(o.maxeval > 0)) throw SnaqError("maxeval must be > 0 (maxeval = " + std::to_string(o.maxeval) + ").");
    if (!(maxequivPLs > 0)) throw SnaqError("maxequivPLs must be > 0 (maxequivPLs = " + std::to_string(maxequivPLs) + ").");
    if (!(0 <= rho && rho <= 1)) throw SnaqError("ρ must be in range [0, 1] (ρ = " + julia_float(rho) + ")");
    if (!(0 < propQuartets && propQuartets <= 1))
        throw SnaqError("propQuartets must be in range (0, 1] (propQuartets = " + julia_float(propQuartets) + ")");
    if (!(0 <= propQuartetsFinal && propQuartetsFinal <= 1))
        throw SnaqError("propQuartetsFinal must be in range [0, 1] (propQuartetsFinal = " +
                        julia_float(propQuartetsFinal) + ")");
    if (!(0 <= o.probQR && o.probQR <= 1)) throw SnaqError("probQR must be in range [0, 1] (probQR = " + julia_float(o.probQR) + ")");
    if (!(0 <= o.probST && o.probST <= 1)) throw SnaqError("probST must be in range [0, 1] (probST = " + julia_float(o.probST) + ")");
    if (o.outgroup != "none") {
        bool found = false;
        for (NodeId l : N.leaf) found |= (N.N(l).name == o.outgroup);
        if (!found) throw SnaqError("No taxa in N have taxa name " + o.outgroup + " (outgroup name)");
    }
    if (!(o.qtolAbs >= 0.0)) throw SnaqError("qtolAbs must be ≥ 0.0 (qtolAbs = " + julia_float(o.qtolAbs) + ")");
    if (!(hmax >= 0)) throw SnaqError("hmax must be ≥ 0 (hmax = " + std::to_string(hmax) + ")");
    if (N.numhybrids > hmax)
        throw SnaqError("The starting network has " + std::to_string(N.numhybrids) +
                        " reticulations, more than hmax = " + std::to_string(hmax) +
                        ". Use `multisearch` or `snaq!` to remove reticulations automatically.");

    const double starttime = wall_time();
    if (!filename.empty()) std::ofstream(filename + ".log", std::ios::trunc);
    if (!filename.empty())
        logmessage(filename, "BEGIN: search with seed " + std::to_string(seed) + " at " + current_time() +
                                 "\n       starting topology: " + writenewick(N, true));

    Rng rng(static_cast<uint64_t>(seed));
    N = N.deepcopynetwork();
    for (NodeId n : N.node)
        if (!N.N(n).leaf && !N.N(n).hybrid) N.N(n).name.clear();
    semidirectnetwork(N);
    if (!o.restrictions(N)) throw SnaqError("N does not meet restrictions IMMEDIATELY");

    if (rng.rand() < o.probST) {
        bool found_different_net = false;
        for (int j = 1; j <= 10000; ++j) {
            try {
                MoveParams p = samplerNNIparameters(N, 1, rng);
                if (!p.valid) throw SnaqError("MethodError: no NNI parameters");
                performrNNI(N, p.p[0].id, p.p[1].id, p.p[2].id, p.p[3].id, 1, rng);
                if (o.restrictions(N)) found_different_net = true;
            } catch (const std::exception&) {
            }
            // QUIRK: as in SNaQ.jl's `finally` block, a failed attempt is not
            // reverted: N is re-copied (with the attempted move) and retried.
            if (!found_different_net) {
                N = N.deepcopynetwork();
                semidirectnetwork(N);
            } else {
                break;
            }
        }
        if (!found_different_net)
            warn("Failed to adjust input network via NNI moves in a manner that met the provided restrictions "
                 "(probST). Using the provided network instead.");
    }

    const long long nqall = nchoose4(N.numtaxa);
    std::vector<long long> q_idxs;
    if (o.qinfTest) {
        const CFMatrix& full = q.matrix();
        std::vector<char> informative(full.nrow(), 1);
        for (size_t i = 0; i < full.nrow(); ++i) informative[i] = isquartetinformative(full.row(i), o.qtolAbs);
        q_idxs = sampleqindices(nqall, propQuartets, informative, rng);
    } else if (propQuartets == 1.0) {
        q_idxs = allranks(nqall);
    } else {
        q_idxs = sampleqindices(nqall, propQuartets, rng);
    }
    double current_logPL = 0.0;
    std::vector<double> CFdeltas;
    long long unchanged_iters = 0;
    CFMatrix qsub = q.rows(q_idxs, nthreads);
    // The search works on quartet classes (classes.hpp): the equations of the
    // classes a move leaves untouched are taken from the current network.
    const QuartetList qlist = QuartetList::fromranks(N.numtaxa, q_idxs);
    const std::vector<double> wsub = o.geneWeights ? q.weights(q_idxs) : std::vector<double>();
    q_idxs = std::vector<long long>();
    ClassOptions copt;
    copt.perquartet = o.probQR > 0.0;
    copt.weights = wsub.empty() ? nullptr : wsub.data();
    copt.qlogq = obsqlogq(qsub, copt.weights);
    copt.nthreads = nthreads;
    fixnegativeedges(N);
    checkgammas(N);
    N.compact_pids();
    QuartetClasses N_cls = buildclasses(N, qlist, qsub, copt);
    // Classes of `net` (a copy of N after a move, or N itself), reusing N's equations.
    auto classesof = [&](Network& net, bool defer = false) {
        ClassOptions po = copt;
        po.reuse = &N_cls;
        po.defer = defer;
        return buildclasses(net, qlist, qsub, po);
    };

    if (o.updateBL) {
        FitInfo fi;
        fitclasses(N, N_cls, rho, search_fit_options(o, -1), nthreads, &fi);
        if (!o.restrictions(N)) throw SnaqError("N does not meet restrictions after preopt");
        current_logPL = N.fscore;
        if (fi.flipped) N_cls = classesof(N);
    } else {
        current_logPL = scoreclasses(N, N_cls, rho, nthreads);
    }

    // Whether a score would be accepted over the current network's.
    // QUIRK: the relative criterion divides by the current score, which is
    // negative; it is undefined (NaN) if the current score is 0.
    auto improves = [&](double score) {
        return score - current_logPL > o.liktolAbs && (current_logPL - score) / current_logPL > o.liktolRel;
    };
    // Topologies this run has evaluated (o.dedup): the score of their fit, or
    // an upper bound of it. An evaluation is final if its fit converged, if it
    // was a bound, or if it was made from the current network (fitting again
    // would repeat the same computation). `epoch` counts accepted moves.
    struct Evaluated {
        double score;
        bool converged;
        long long epoch;
    };
    std::unordered_map<uint64_t, Evaluated> evaluated;
    long long epoch = 0, dupskips = 0;
    auto record = [&](uint64_t sig, double score, bool converged) {
        auto it = evaluated.find(sig);
        if (it == evaluated.end()) {
            evaluated.emplace(sig, Evaluated{score, converged, epoch});
        } else {
            it->second.score = std::max(it->second.score, score);
            it->second.converged = it->second.converged || converged;
            it->second.epoch = epoch;
        }
        if (shared && converged) {
            std::lock_guard<std::mutex> g(shared->m);
            auto sit = shared->score.find(sig);
            if (sit == shared->score.end())
                shared->score.emplace(sig, score);
            else
                sit->second = std::max(sit->second, score);
        }
    };
    if (o.dedup) {
        bool discrete = false;
        const uint64_t sig = topologysignature(N, &discrete);
        if (discrete) evaluated.emplace(sig, Evaluated{current_logPL, true, epoch});
    }

    std::vector<AttemptedMove> moves_attempted;
    MoveStats stats;
    // Hybrid additions (guideHybrids, minAddHybrid)
    const long long minaddhybrid = o.minAddHybrid >= 0 ? o.minAddHybrid : maxequivPLs;
    long long addhybrid_evaluated = 0;  // since the last accepted move
    long long sweep_skipped = 0;        // proposals in a row that were not evaluated, while only adding hybrids
    bool hybridsweep = false;           // only hybrid additions are proposed
    std::vector<HybridCandidate> candidates;
    bool candidates_current = false;
    logtext(logfile, "Entering main loop with -logPL = " + julia_float(current_logPL));

    for (long long j = 2; j <= o.maxeval; ++j) {
        if (j % 100 == 0) logmoves(logfile, stats);
        if (o.verbose) {
            std::lock_guard<std::mutex> g(g_stdout_mutex);
            std::cout << "\rIteration " << j << "/" << o.maxeval << " - in a row=" << unchanged_iters << "/"
                      << maxequivPLs << "              " << std::flush;
        }
        Network Nprime = N.deepcopynetwork();
        ProposalGuide guide;
        guide.forceaddhybrid = hybridsweep;
        if (o.guideHybrids && N.numhybrids < hmax) {
            if (!candidates_current) {
                candidates = hybridcandidates(N, N_cls, rho);
                candidates_current = true;
            }
            guide.candidates = &candidates;
        }
        std::pair<Move, MoveParams> proposal =
            generatemoveproposal(Nprime, &N_cls, moves_attempted, hmax, o.probQR, qsub, CFdeltas, rng, rho, &guide);
        const Move prop_move = proposal.first;
        const MoveParams prop_params = proposal.second;
        if (prop_move == Move::none) {
            logtext(logfile, "Iteration " + std::to_string(j) + ": no untried move proposals remain, stopping.");
            break;
        }
        applymove(Nprime, prop_move, prop_params, rng);
        const std::string mname = move_name(prop_move);
        stats.proposed[mname] += 1;
        stats.accepted[mname] += 0;
        stats.dPL[mname];

        bool accepted = false;
        bool duplicate = false;  // the proposal is a topology already evaluated (o.dedup)
        auto prefix = [&]() {
            return "Iteration " + std::to_string(j) + " (N.h=" + std::to_string(N.numhybrids) +
                   "), in a row = " + std::to_string(unchanged_iters) + "/" + std::to_string(maxequivPLs);
        };
        // body of the iteration; `false` means `continue`
        auto body = [&]() -> bool {
            removedegree2nodes(Nprime);
            while (shrink3cycles(Nprime)) {
            }
            while (shrink2cycles(Nprime)) {
            }
            if (o.outgroup != "none") {
                try {
                    rootatnode_name(Nprime, o.outgroup);
                } catch (const RootMismatch&) {
                    logtext(logfile, prefix() + " REJECTED " + mname + " (cannot reroot at outgroup)");
                    return false;
                }
            }
            semidirectnetwork(Nprime);
            if (!o.restrictions(Nprime)) {
                logtext(logfile, prefix() + " REJECTED " + mname + " (restrictions not met)");
                return false;
            }
            uint64_t sig = 0;
            bool havesig = false;
            if (o.dedup) {
                sig = topologysignature(Nprime, &havesig);
                if (havesig) {
                    auto it = evaluated.find(sig);
                    if (it != evaluated.end() && (it->second.converged || it->second.epoch == epoch) &&
                        !improves(it->second.score)) {
                        logtext(logfile, prefix() + " REJECTED " + mname + " (topology already evaluated)");
                        if (o.dedupCounts) {
                            unchanged_iters += 1;
                            return true;
                        }
                        // not a new proposal: draw another one, as for moves already attempted
                        ++dupskips;
                        duplicate = true;
                        return false;
                    }
                    if (shared) {
                        bool known = false;
                        double other = 0.0;
                        {
                            std::lock_guard<std::mutex> g(shared->m);
                            auto sit = shared->score.find(sig);
                            if (sit != shared->score.end()) {
                                known = true;
                                other = sit->second;
                            }
                        }
                        if (known && !improves(other)) {
                            unchanged_iters += 1;
                            logtext(logfile, prefix() + " REJECTED " + mname + " (" + round3(other) + " < " +
                                                 round3(current_logPL) + ", fitted by another run)");
                            evaluated.emplace(sig, Evaluated{other, true, epoch});
                            return true;
                        }
                    }
                }
            }
            dupskips = 0;
            // (SNaQ.jl updates the equations in place after rNNI moves and rebuilds
            // them otherwise; here every move reuses the equations of the classes
            // it leaves untouched.)
            QuartetClasses Nprime_cls = classesof(Nprime, o.boundPrune);
            if (o.boundPrune) {
                // no fit of Nprime can score more than `bound`
                const double bound = classupperbound(Nprime_cls);
                if (!improves(bound)) {
                    unchanged_iters += 1;
                    logtext(logfile, prefix() + " REJECTED " + mname + " (bound " + round3(bound) + " < " +
                                         round3(current_logPL) + ")");
                    if (havesig) record(sig, bound, true);
                    return true;
                }
                finishclasses(Nprime, Nprime_cls, nthreads);
            }
            FitInfo fi;
            double Nprime_logPL =
                fitclasses(Nprime, Nprime_cls, rho, search_fit_options(o, opt_maxeval), nthreads, &fi);
            if (Nprime_logPL == -std::numeric_limits<double>::infinity())
                throw SnaqError("Nprime_logPL is -Inf?? newick: " + writenewick(Nprime, true));
            if (std::isnan(Nprime_logPL))
                throw SnaqError("Nprime_logPL = NaN\n" + writenewick(Nprime, true));
            if (improves(Nprime_logPL)) {
                N = std::move(Nprime);
                N_cls = std::move(Nprime_cls);
                // equations refer to a hybrid's minor edge: rebuild those of the
                // hybrids whose major and minor edges were exchanged by the fit
                if (fi.flipped) N_cls = classesof(N);
                CFdeltas.clear();
                stats.accepted[mname] += 1;
                stats.dPL[mname].push_back(Nprime_logPL - current_logPL);
                current_logPL = Nprime_logPL;
                logtext(logfile, "Iteration " + std::to_string(j) + " (N.h=" + std::to_string(N.numhybrids) +
                                     "), in a row = " + std::to_string(unchanged_iters) + "/" +
                                     std::to_string(maxequivPLs) + " ACCEPTED " + mname +
                                     ", new -logPL=" + round6(current_logPL));
                unchanged_iters = 0;
                moves_attempted.clear();
                accepted = true;
                ++epoch;
                candidates_current = false;
                addhybrid_evaluated = 0;
                hybridsweep = false;
                // a topology that was the current network is never fitted again
                if (havesig) record(sig, current_logPL, true);
            } else {
                unchanged_iters += 1;
                logtext(logfile, prefix() + " REJECTED " + mname + " (" + round3(Nprime_logPL) + " < " +
                                     round3(current_logPL) + ")");
                if (havesig) record(sig, Nprime_logPL, fi.converged);
            }
            return true;
        };
        bool completed = body();
        // (duplicates are not recorded: the list would grow without bound, and
        // the search stops anyway after 1000 duplicates in a row)
        if (!accepted && !duplicate) moves_attempted.push_back(snapshot_attempt(Nprime, prop_move, prop_params));
        if (completed && !accepted && prop_move == Move::addhybrid) ++addhybrid_evaluated;
        if (hybridsweep) sweep_skipped = completed ? 0 : sweep_skipped + 1;
        if (unchanged_iters > maxequivPLs && (completed || hybridsweep)) {
            // The run would stop here. With room for more reticulations, first
            // make sure that enough hybrid additions were evaluated.
            if (N.numhybrids < hmax && addhybrid_evaluated < minaddhybrid && sweep_skipped < 1000) {
                if (!hybridsweep)
                    logtext(logfile, "Iteration " + std::to_string(j) + ": " + std::to_string(addhybrid_evaluated) +
                                         " hybrid additions evaluated since the last accepted move; proposing only "
                                         "hybrid additions until " + std::to_string(minaddhybrid) + " are.");
                hybridsweep = true;
            } else {
                break;
            }
        }
        if (dupskips >= 1000) {
            logtext(logfile, "Iteration " + std::to_string(j) +
                                 ": 1000 proposals in a row were topologies already evaluated, stopping.");
            break;
        }
    }
    if (o.verbose) {
        std::lock_guard<std::mutex> g(g_stdout_mutex);
        std::cout << "\n";
    }
    N.fscore = current_logPL;

    FitOptions finalfit;
    finalfit.method = fit_method(o);
    if (propQuartetsFinal == 1.0 && propQuartets != 1.0) {
        logmessage(filename, "Re-optimizing branch lengths with ALL quartets.");
        // QUIRK: SNaQ.jl calls fitnumericalparameters!(N, q) here, which uses ρ = 0
        // regardless of the ρ of the search.
        const CFMatrix& full = q.matrix();
        const std::vector<double> fw = o.geneWeights ? q.weights(allranks(nqall)) : std::vector<double>();
        ClassOptions fo;
        fo.weights = fw.empty() ? nullptr : fw.data();
        fo.qlogq = obsqlogq(full, fo.weights);
        fo.nthreads = nthreads;
        fixnegativeedges(N);
        checkgammas(N);
        QuartetClasses fc = buildclasses(N, QuartetList::allquartets(N.numtaxa), full, fo);
        N.fscore = fitclasses(N, fc, 0.0, finalfit, nthreads);
        logmessage(filename,
                   "END propQuartets<1.0 post-search parameter optimization: found minimizer topology with SNaQ "
                   "score=" + round5(N.fscore));
    } else if (0 < propQuartetsFinal && propQuartetsFinal < 1) {
        std::vector<long long> final_idxs = sampleqindices(nchoose4(N.numtaxa), propQuartetsFinal, rng);
        logmessage(filename, "Re-optimizing branch lengths with " + std::to_string(final_idxs.size()) + " of " +
                                 std::to_string(nchoose4(N.numtaxa)) + " quartets (propQuartetsFinal = " +
                                 julia_float(propQuartetsFinal) + ").");
        const CFMatrix fq = q.rows(final_idxs, nthreads);
        const std::vector<double> fw = o.geneWeights ? q.weights(final_idxs) : std::vector<double>();
        ClassOptions fo;
        fo.weights = fw.empty() ? nullptr : fw.data();
        fo.qlogq = obsqlogq(fq, fo.weights);
        fo.nthreads = nthreads;
        fixnegativeedges(N);
        checkgammas(N);
        QuartetClasses fc = buildclasses(N, QuartetList::fromranks(N.numtaxa, final_idxs), fq, fo);
        N.fscore = fitclasses(N, fc, rho, finalfit, nthreads);
        logmessage(filename,
                   "END propQuartetsFinal<1.0 post-search parameter optimization: found minimizer topology with "
                   "SNaQ score=" + round5(N.fscore));
    }

    for (NodeId n : N.node)
        if (!N.N(n).leaf && !N.N(n).hybrid) N.N(n).name.clear();
    for (size_t iH = 0; iH < N.hybrid.size(); ++iH) N.N(N.hybrid[iH]).name = "H" + std::to_string(iH + 1);

    logtext(logfile, "Search complete at " + current_time() + ".\n\n");
    logmoves(logfile, stats);
    if (!filename.empty()) {
        logmessage(filename, "END: search with seed " + std::to_string(seed) + " after " +
                                 time_elapsed(wall_time() - starttime) + ". SNaQ score = " + julia_float(N.fscore));
        logmessage(filename, writenewick(N));
    }
    return N;
}

// ------------------------------------------------------------------------
// multisearch
// ------------------------------------------------------------------------

static bool approx_eq(double a, double b) {
    return std::abs(a - b) <= 1.4901161193847656e-8 * std::max(std::abs(a), std::abs(b));
}

static std::vector<Network> verifystartingtopologies(const std::vector<Network>& N, bool single_input,
                                                     const std::string& outgroup, const Restrictions& restrictions,
                                                     int hmax) {
    std::vector<Network> Ns;
    for (const Network& n : N) Ns.push_back(n.deepcopynetwork());
    for (size_t j = 0; j < Ns.size(); ++j) {
        Network& net = Ns[j];
        const std::string num = std::to_string(j + 1);
        bool poly = false;
        for (NodeId n : net.node) poly |= net.N(n).edge.size() > 3;
        if (poly) {
            warn("Input network #" + num +
                 " has a polytomy. SNaQ only infers binary networks, so this will be automatically resolved before "
                 "inference.");
            // QUIRK: SNaQ.jl bounds this loop with `N.numnodes`, where N is the
            // user's input: this errors when several starting topologies are given.
            if (!single_input) throw SnaqError("type Array has no field numnodes");
            int iters = 0;
            while (true) {
                ++iters;
                if (iters > N[0].numnodes)
                    throw SnaqError("Got stuck in an infinite loop while resolving polytomies. Please report this bug "
                                    "with your input tree(s) on GitHub.");
                NodeId multi = NONE;
                for (NodeId n : net.node)
                    if (net.N(n).edge.size() > 3) {
                        multi = n;
                        break;
                    }
                if (multi == NONE) break;
                resolvetreepolytomy(net, multi);
            }
        }
        semidirectnetwork(net);
        for (EdgeId e : net.edge)
            if (net.E(e).length == -1.0) net.E(e).length = 0.0;
        for (NodeId h : net.hybrid) {
            double gM = net.E(net.getparentedge(h)).gamma, gm = net.E(net.getparentedgeminor(h)).gamma;
            if (1 >= gM && gM >= 0 && 1 >= gm && gm >= 0 && approx_eq(gM + gm, 1.0)) continue;
            net.E(net.getparentedge(h)).gamma = 0.5;
            net.E(net.getparentedgeminor(h)).gamma = 0.5;
        }
        if (net.numhybrids > hmax) {
            warn("Starting topology #" + num + " has " + std::to_string(net.numhybrids) +
                 " reticulations but hmax = " + std::to_string(hmax) + ": removing " +
                 std::to_string(net.numhybrids - hmax) +
                 " of them, the reticulations with the smallest minor inheritance probability γ first.");
            removeexcessreticulations(net, hmax);
        }
        if (!restrictions(net))
            throw ArgumentError("Starting topology #" + num + " does not meet provided restrictions (" +
                                restrictions.describe() + ").");
        if (outgroup == "none") continue;
        bool has = false;
        for (NodeId l : net.leaf) has |= net.N(l).name == outgroup;
        if (!has)
            throw ArgumentError("Starting topology #" + num + " does not contain the supplied outgroup (" + outgroup +
                                ").");
        try {
            rootatnode_name(net, outgroup);
        } catch (const RootMismatch&) {
            throw ArgumentError("Starting topology #" + num +
                                " contains the outgroup but cannot be rooted at the outgroup.");
        }
    }
    return Ns;
}

static std::string abspath(const std::string& p) {
    return fs::absolute(fs::path(p)).lexically_normal().string();
}

long long rescoreonsharedquartets(std::vector<Network>& nets, const CFSource& q, double propQuartetsFinal,
                                  long long seed, double rho, int nthreads, bool geneweights) {
    Rng rng(static_cast<uint64_t>(seed));
    std::vector<long long> idxs = sampleqindices(q.nrows(), propQuartetsFinal, rng);
    CFMatrix qsub = q.rows(idxs, nthreads);
    const std::vector<double> w = geneweights ? q.weights(idxs) : std::vector<double>();
    ClassOptions co;
    co.weights = w.empty() ? nullptr : w.data();
    co.qlogq = obsqlogq(qsub, co.weights);
    co.nthreads = nthreads;
    for (Network& net : nets) {
        fixnegativeedges(net);
        checkgammas(net);
        QuartetClasses cls = buildclasses(net, QuartetList::fromranks(net.numtaxa, idxs), qsub, co);
        net.fscore = scoreclasses(net, cls, rho, nthreads);
    }
    return static_cast<long long>(idxs.size());
}

std::vector<Network> multisearch(std::vector<Network> N, bool single_input, const CFSource& q, int hmax,
                                 const SnaqOptions& o, double propQuartetsFinal) {
    const int runs = o.runs;
    if (!(runs > 0)) throw SnaqError("runs must be > 0 (runs = " + std::to_string(runs) + ").");
    if (!single_input && N.size() != 1 && static_cast<int>(N.size()) != runs)
        throw SnaqError("If N is a vector, it must be length 1 or length equal to runs (length(N) = " +
                        std::to_string(N.size()) + ", runs = " + std::to_string(runs) + ")");

    if (!q.is_lazy()) {
        // taxa of the data must all be in each input; input taxa not in the data are deleted
        const std::vector<std::string>& dcf_names = q.taxa();
        for (Network& T : N) {
            auto tl = T.tiplabels();
            std::set<std::string> tset(tl.begin(), tl.end());
            for (const auto& t : dcf_names)
                if (!tset.count(t)) throw SnaqError("Taxa in DataCF does not match taxa in input.");
            std::set<std::string> dset(dcf_names.begin(), dcf_names.end());
            std::vector<std::string> removetaxa;
            for (const auto& t : tl)
                if (!dset.count(t)) removetaxa.push_back(t);
            if (!removetaxa.empty()) {
                std::string list;
                for (size_t k = 0; k < removetaxa.size(); ++k) list += (k ? ", " : "") + ("\"" + removetaxa[k] + "\"");
                warn("The following are in the inputs but not the DataCF, these taxa will be deleted: Any[" + list + "]");
                for (const auto& t : removetaxa) deleteleaf_name(T, t, DeleteOpts{});
            }
        }
        for (Network& T : N) {
            auto tl = T.tiplabels();
            std::set<std::string> a(tl.begin(), tl.end()), b(dcf_names.begin(), dcf_names.end());
            if (a != b) throw SnaqError("Taxa in DataCF does not match taxa in input.");
        }
    } else {
        long long nrow = nchoose4(N[0].numtaxa);
        if (q.nrows() != nrow)
            throw SnaqError("Input CF matrix should have size (" + std::to_string(nrow) + ", 3), has size (" +
                            std::to_string(q.nrows()) + ", 3) instead.");
    }

    if (!(hmax >= 0)) throw SnaqError("hmax must be ≥ 0 (hmax = " + std::to_string(hmax) + ").");
    std::vector<Network> Ns = verifystartingtopologies(N, single_input, o.outgroup, o.restrictions, hmax);

    Rng seedrng(static_cast<uint64_t>(o.seed));
    std::vector<long long> run_seeds(static_cast<size_t>(runs));
    for (auto& s : run_seeds) s = std::llabs(seedrng.rand_int64() % 100000);

    const std::string& filename = o.filename;
    const int nthreads_total = std::max(1, o.threads);
    std::string msg = "Beginning network optimization using SNaQ.jl across " + std::to_string(runs) +
                      " runs with the following parameters:\n    hmax = " + std::to_string(hmax) +
                      ",\n    seed = " + std::to_string(o.seed) + ",\n    outgroup = " + o.outgroup +
                      ",\n    restrictions = " + o.restrictions.describe() + "\nRoot name for log files: " + filename +
                      " (absolute path " + abspath(filename) + ")\nCurrently utilizing 1 processor and " +
                      std::to_string(nthreads_total) + " thread" + (nthreads_total > 1 ? "s" : "") + ".\n";
    if (o.verbose) std::cout << msg << std::endl;
    std::vector<std::string> filenames(static_cast<size_t>(runs));
    if (!filename.empty()) {
        std::ofstream(filename + ".log", std::ios::trunc);
        logmessage(filename, msg);
        std::string runs_path = filename + "_runs/";
        fs::create_directories(runs_path);
        for (int j = 0; j < runs; ++j) filenames[static_cast<size_t>(j)] = runs_path + "run" + std::to_string(j + 1);
    }

    // Runs are independent: spread them over the threads; leftover threads
    // parallelize the quartet loops inside each run.
    const double starttime = wall_time();
    std::vector<Network> all_nets(static_cast<size_t>(runs));
    const int concurrent = std::min(runs, nthreads_total);
    const int inner = std::max(1, nthreads_total / concurrent);
    // With propQuartets = 1 every run fits to the same quartets, so a topology
    // fitted by one run need not be fitted again by another.
    SearchShared sharedscores;
    SearchShared* sharedptr = (o.dedup && o.propQuartets == 1.0) ? &sharedscores : nullptr;
    std::atomic<int> next{0};
    std::exception_ptr err = nullptr;
    std::mutex errm;
    auto worker = [&]() {
        while (true) {
            int j = next.fetch_add(1);
            if (j >= runs) return;
            {
                std::lock_guard<std::mutex> g(errm);
                if (err) return;
            }
            try {
                const Network& start = Ns.size() == 1 ? Ns[0] : Ns[static_cast<size_t>(j)];
                std::string logfile = o.logprefix.empty() ? "" : o.logprefix + std::to_string(j + 1);
                all_nets[static_cast<size_t>(j)] =
                    search(start.clone(), q, hmax, o, run_seeds[static_cast<size_t>(j)], propQuartetsFinal, logfile,
                           filenames[static_cast<size_t>(j)], inner, sharedptr);
            } catch (...) {
                std::lock_guard<std::mutex> g(errm);
                if (!err) err = std::current_exception();
            }
        }
    };
    std::vector<std::thread> pool;
    for (int t = 1; t < concurrent; ++t) pool.emplace_back(worker);
    worker();
    for (auto& t : pool) t.join();
    if (err) std::rethrow_exception(err);
    const std::string elapsed = time_elapsed(wall_time() - starttime);

    if (0 < propQuartetsFinal && propQuartetsFinal < 1) {
        long long nused =
            rescoreonsharedquartets(all_nets, q, propQuartetsFinal, o.seed, o.rho, nthreads_total, o.geneWeights);
        if (!filename.empty())
            logmessage(filename, "Re-scored the networks of all " + std::to_string(runs) + " runs on the same " +
                                     std::to_string(nused) + " quartets (propQuartetsFinal = " +
                                     julia_float(propQuartetsFinal) + ").");
    }

    std::vector<size_t> sort_idx(all_nets.size());
    std::iota(sort_idx.begin(), sort_idx.end(), size_t{0});
    std::stable_sort(sort_idx.begin(), sort_idx.end(),
                     [&](size_t a, size_t b) { return all_nets[a].fscore > all_nets[b].fscore; });
    Network& bestnet = all_nets[sort_idx[0]];

    if (!filename.empty()) {
        logmessage(filename, "Finished optimizing topology at " + current_time() + " after " + elapsed +
                                 ".\nOptimal network: " + writenewick(bestnet, true) + "\nOptimal SNaQscore: " +
                                 julia_float(bestnet.fscore) + "\nTo view all " + std::to_string(runs) +
                                 " inferred networks and their associated SNaQscore scores, see " + filename +
                                 ".out (" + abspath(filename + ".out") + ")");
        {
            std::ofstream f(filename + ".out", std::ios::trunc);
            f << writenewick(bestnet) << " SNaQ score = " << julia_float(bestnet.fscore) << "\nElapsed time: " << elapsed
              << ", " << runs << " attempted runs\n\n-----------------------------------\n"
              << "List of estimated networks for all runs (sorted by log-pseudolik; the larger, the better):\n";
            for (size_t j : sort_idx)
                f << " " << writenewick(all_nets[j]) << ", with SNaQ score " << julia_float(all_nets[j].fscore) << "\n";
            f << "-----------------------------------\n";
        }
        {
            std::ofstream f(filename + ".networks", std::ios::trunc);
            for (size_t k = 0; k < sort_idx.size(); ++k) {
                size_t i = sort_idx[k];
                f << writenewick(all_nets[i]) << ", with SNaQ score " << julia_float(all_nets[i].fscore);
                if (k == 0) f << " (best network found, remaining sorted by log-pseudolik; the larger, the better)";
                f << "\n";
            }
        }
    }
    for (Network& n : all_nets)
        for (NodeId l : n.leaf) n.E(n.getparentedge(l)).length = -1;

    std::vector<Network> sorted;
    sorted.reserve(all_nets.size());
    for (size_t i : sort_idx) sorted.push_back(std::move(all_nets[i]));
    return sorted;
}

}  // namespace csnaq
