// Validation tool: compares the quartet-class likelihood (classes.hpp) with
// the per-quartet likelihood on every network of a Newick file: SNaQ score
// and gradient at the network's own parameter values, for ρ = 0 and ρ = 0.5,
// against synthetic observed CFs.
//
// It also checks expected CFs at ρ = 0 of some quartets of each network
// against an independent brute-force calculation.
//
// usage: csnaq_classcheck <networks file> [max quartets per network] [first line] [last line]
//                         [brute-force quartets per network]
// prints one line per network and a summary; exits with 1 if a class's expected CFs differ from
// the brute-force values, or if the equations of the two builders (porteqn.hpp) differ in value. (Scores can differ between the two likelihoods: see the summary.)
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "bruteforce_cf.hpp"
#include "classes.hpp"
#include "equations.hpp"
#include "likelihood.hpp"
#include "newick.hpp"
#include "porteqn.hpp"
#include "quartets.hpp"
#include "snaq_network.hpp"

using namespace csnaq;

static void synthetic_obs(long long r, double* row) {
    double a = 1 + ((r * 7) % 11), b = 1 + ((r * 13 + 5) % 7), c = 1 + ((r * 3 + 1) % 5);
    if (r % 4 == 0) c = 0;
    if (r % 9 == 0) a = 0;
    double s = a + b + c;
    row[0] = a / s;
    row[1] = b / s;
    row[2] = c / s;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: csnaq_classcheck <networks file> [max quartets] [first line] [last line]\n";
        return 2;
    }
    const long long cap = argc > 2 ? std::stoll(argv[2]) : 20000;
    const long first = argc > 3 ? std::stol(argv[3]) : 1;
    const long last = argc > 4 ? std::stol(argv[4]) : 1000000000;
    std::ifstream in(argv[1]);
    std::string s;
    long line = 0, checked = 0, bad = 0, errors = 0, tot_ochecked = 0, tot_obad_old = 0, tot_obad_new = 0;
    const long noracle = argc > 5 ? std::stol(argv[5]) : 200;
    double worst = 0.0, worstref = 0.0;
    long refbad = 0;
    size_t tot_size_new = 0, tot_size_ref = 0;
    while (std::getline(in, s)) {
        ++line;
        if (line < first || line > last || s.find('(') == std::string::npos) continue;
        std::string olderr, newerr;
        long obad_old = 0, obad_new = 0, ochecked = 0;
        double sold[2] = {0, 0}, snew[2] = {0, 0}, gdiff = 0.0, refdiff = 0.0;
        size_t size_new = 0, size_ref = 0;
        size_t nclass = 0, nq = 0;
        int ntaxa = 0, nhyb = 0;
        std::vector<long long> ranks;
        CFMatrix obs;
        try {
            Network net = readnewick_string(s).deepcopynetwork();
            semidirectnetwork(net);
            ntaxa = net.numtaxa;
            nhyb = net.numhybrids;
            const long long nqall = nchoose4(net.numtaxa);
            if (cap <= 0 || nqall <= cap) {
                for (long long i = 0; i < nqall; ++i) ranks.push_back(i);
            } else {
                for (long long i = 0; i < cap; ++i) ranks.push_back((i * nqall) / cap);
            }
            nq = ranks.size();
            obs.v.resize(3 * nq);
            for (size_t j = 0; j < nq; ++j) synthetic_obs(ranks[j], &obs.v[3 * j]);
            std::vector<double> gold[2], gnew[2];
            std::vector<int> pidof;
            try {
                Network a = net.deepcopynetwork();
                EquationSet es = findquartetequations(a, ranks, 1);
                for (int k = 0; k < 2; ++k)
                    sold[k] = computelossandgradient(es.eqns, es.info.params, gold[k], obs, rhotoalpha(k ? 0.5 : 0.0), 1);
                for (const ParamObj& o : es.info.idx_obj) pidof.push_back(o.is_hybrid_node ? a.N(o.id).pid : a.E(o.id).pid);
            } catch (const std::exception& e) {
                olderr = e.what();
            }
            try {
                Network b = net.deepcopynetwork();
                fixnegativeedges(b, false);
                ClassOptions co;
                co.qlogq = obsqlogq(obs);
                QuartetClasses cls = buildclasses(b, QuartetList::fromranks(b.numtaxa, ranks), obs, co);
                nclass = cls.size();
                std::vector<double> params = pidparams(b, cls.npid), g(static_cast<size_t>(cls.npid));
                LossWorkspace ws;
                for (int k = 0; k < 2; ++k) {
                    snew[k] = evalclasses(cls.batch, params.data(), g.data(), rhotoalpha(k ? 0.5 : 0.0), ws, 1);
                    gnew[k] = g;
                }
                // the same classes with the reference equation builder: same score and gradient
                {
                    Network r = net.deepcopynetwork();
                    fixnegativeedges(r, false);
                    ClassOptions ro = co;
                    ro.referenceequations = true;
                    QuartetClasses ref = buildclasses(r, QuartetList::fromranks(r.numtaxa, ranks), obs, ro);
                    if (ref.size() != cls.size() || ref.npid != cls.npid) throw SnaqError("reference classes differ");
                    std::vector<double> gr(static_cast<size_t>(ref.npid));
                    for (int k = 0; k < 2; ++k) {
                        const double sr = evalclasses(ref.batch, params.data(), gr.data(), rhotoalpha(k ? 0.5 : 0.0), ws, 1);
                        refdiff = std::max(refdiff, std::abs(sr - snew[k]) / std::max(1.0, std::abs(sr)));
                        for (size_t j = 0; j < gr.size(); ++j)
                            refdiff = std::max(refdiff, std::abs(gr[j] - gnew[k][j]) / std::max(1.0, std::abs(gr[j])));
                    }
                    for (size_t ci = 0; ci < cls.size(); ++ci) {
                        size_new += equationsize(cls.eqn[ci]->eqn);
                        size_ref += equationsize(ref.eqn[ci]->eqn);
                    }
                }
                // expected CFs at ρ = 0 against the brute-force oracle, for some quartets
                {
                    Network a = net.deepcopynetwork();
                    EquationSet es = findquartetequations(a, ranks, 1);
                    ClassOptions pq = co;
                    pq.perquartet = true;
                    Network b2 = net.deepcopynetwork();
                    fixnegativeedges(b2, false);
                    QuartetClasses cq = buildclasses(b2, QuartetList::fromranks(b2.numtaxa, ranks), obs, pq);
                    std::vector<double> pp = pidparams(b2, cq.npid);
                    BruteForceCF oracle(b2);
                    std::vector<std::pair<std::string, NodeId>> lv;
                    for (NodeId l : b2.leaf) lv.emplace_back(b2.N(l).name, l);
                    std::sort(lv.begin(), lv.end());
                    const size_t step = std::max<size_t>(1, nq / static_cast<size_t>(noracle));
                    for (size_t j = 0; j < nq; j += step) {
                        auto t = unrank4taxa(b2.numtaxa, ranks[j]);
                        double want[3];
                        if (!oracle.ecf({lv[static_cast<size_t>(t[0])].second, lv[static_cast<size_t>(t[1])].second,
                                         lv[static_cast<size_t>(t[2])].second, lv[static_cast<size_t>(t[3])].second},
                                        want))
                            continue;
                        ++ochecked;
                        auto eo = computeexpectedCF(*es.eqns[j], es.info.params, 0.0);
                        auto ec = classexpectedCF(*cq.eqn[cq.qclass[j]], pp.data(), cq.npid, rhotoalpha(0.0));
                        double en[3];
                        for (int k = 0; k < 3; ++k) en[(cq.qcols[j] >> (2 * k)) & 3] = ec[static_cast<size_t>(k)];
                        if (std::abs(eo[0] - want[0]) + std::abs(eo[1] - want[1]) > 1e-7) ++obad_old;
                        if (std::abs(en[0] - want[0]) + std::abs(en[1] - want[1]) > 1e-7) ++obad_new;
                    }
                }
                // equations found again in the classes of a copy must give the same result
                Network c = b.deepcopynetwork();
                co.reuse = &cls;
                QuartetClasses cls2 = buildclasses(c, QuartetList::fromranks(c.numtaxa, ranks), obs, co);
                if (cls2.nbuilt != 0) newerr = "equations were rebuilt for an unchanged copy";
                std::vector<double> g2(static_cast<size_t>(cls.npid));
                double s2 = evalclasses(cls2.batch, params.data(), g2.data(), rhotoalpha(0.5), ws, 1);
                if (std::abs(s2 - snew[1]) > 1e-10 * std::max(1.0, std::abs(s2))) newerr = "reuse changed the score";
            } catch (const std::exception& e) {
                newerr = e.what();
            }
            if (olderr.empty() && newerr.empty())
                for (int k = 0; k < 2; ++k)
                    for (size_t j = 0; j < pidof.size(); ++j)
                        gdiff = std::max(gdiff, std::abs(gold[k][j] - gnew[k][static_cast<size_t>(pidof[j])]) /
                                                    std::max(1.0, std::abs(gold[k][j])));
        } catch (const std::exception& e) {
            olderr = newerr = e.what();
        }
        tot_ochecked += ochecked;
        tot_obad_old += obad_old;
        tot_obad_new += obad_new;
        if (!olderr.empty() || !newerr.empty()) {
            const bool both = !olderr.empty() && !newerr.empty();
            if (both) {
                ++errors;
            } else {
                ++bad;
                std::cout << line << " MISMATCH old error: [" << olderr << "] new error: [" << newerr << "]\n";
            }
            continue;
        }
        double d = gdiff;
        for (int k = 0; k < 2; ++k) d = std::max(d, std::abs(sold[k] - snew[k]) / std::max(1.0, std::abs(sold[k])));
        if (std::getenv("CLASSCHECK_VERBOSE"))
            std::cout << "  score rho=0: " << sold[0] << " vs " << snew[0] << "; rho=0.5: " << sold[1] << " vs " << snew[1]
                      << "; gradient reldiff " << gdiff << "\n";
        ++checked;
        worst = std::max(worst, d);
        worstref = std::max(worstref, refdiff);
        tot_size_new += size_new;
        tot_size_ref += size_ref;
        if (refdiff > 1e-9) {
            ++refbad;
            std::cout << line << " REFERENCE MISMATCH reldiff=" << refdiff << " ";
        }
        if (d > 1e-8) {
            ++bad;
            std::cout << line << " MISMATCH ";
        } else {
            std::cout << line << " ok ";
        }
        std::cout << "ntaxa=" << ntaxa << " h=" << nhyb << " quartets=" << nq << " classes=" << nclass
                  << " score=" << sold[0] << " reldiff=" << d << " oracle: " << ochecked << " quartets, "
                  << obad_old << " wrong per-quartet, " << obad_new << " wrong classes; equation nodes: " << size_new
                  << " (reference builder: " << size_ref << ")\n";
    }
    std::cout << "checked " << checked << " network(s), " << errors << " with an error in both implementations, " << bad
              << " with a score or gradient mismatch; worst relative difference " << worst << "\n"
              << "brute-force expected CFs (rho = 0) on " << tot_ochecked << " quartets: " << tot_obad_old
              << " wrong with per-quartet equations, " << tot_obad_new << " wrong with quartet classes\n"
              << "equation builder against the reference builder (score and gradient, rho = 0 and 0.5): " << refbad
              << " network(s) differ; worst relative difference " << worstref << "; equation nodes " << tot_size_new
              << " (reference: " << tot_size_ref << ")\n";
    return (tot_obad_new || refbad) ? 1 : 0;
}
