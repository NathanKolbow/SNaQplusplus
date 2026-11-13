# One parameter-optimization benchmark case with SNaQ.jl's
# fitnumericalparameters! (called by benchmarks/bench.py): fit the branch
# lengths and γ of the true network, from fixed starting values (internal
# branch lengths 1, minor γ 0.3), to the observed CFs of the gene trees.
# Mirrors tools/bench_fit.cpp.
#
# usage: julia --project=tools/julia run_fit_case.jl <truth network> <genetrees> <maxeval> <repeats>
#
# Prints "RESULT {json}" with the best (minimum) times over the repeats.
using SNaQ, PhyloNetworks, DataFrames

base = Sys.maxrss()
netfile, gtfile = ARGS[1], ARGS[2]
maxeval, repeats = parse(Int, ARGS[3]), parse(Int, ARGS[4])
net0 = readnewick(readlines(netfile)[1])
trees = readmultinewick(gtfile)
q4, t4 = countquartetsintrees(trees; showprogressbar=false)
d = readtableCF(DataFrame(tablequartetCF(q4, t4)))
q = SNaQ.gatherCFmatrix(d)

function fitonce()
    net = SNaQ.deepcopynetwork(net0)
    SNaQ.semidirectnetwork!(net)
    for e in net.edge
        getchild(e).leaf || (e.length = 1.0)
    end
    for h in net.hybrid
        getparentedge(h).gamma = 0.7
        getparentedgeminor(h).gamma = 0.3
    end
    teq = @elapsed eqns = SNaQ.findquartetequations(net)[1]
    tfit = @elapsed score = SNaQ.fitnumericalparameters!(net, eqns, q, 0.0; maxeval=maxeval)
    return teq, tfit, score, SNaQ.countoptimizationparams(net), length(eqns)
end

fitonce()  # warm-up: compilation
best_eq = Inf; best_fit = Inf; score = NaN; np = 0; nq = 0
for r in 1:repeats
    teq, tfit, s, n, m = fitonce()
    global best_eq = min(best_eq, teq); global best_fit = min(best_fit, tfit)
    global score = s; global np = n; global nq = m
end
println("RESULT {\"fit_s\": $best_fit, \"equations_s\": $best_eq, \"score\": $score, \"nparams\": $np, ",
        "\"nquartets\": $nq, \"peak_rss_mb\": $(Sys.maxrss() / 2^20), \"base_rss_mb\": $(base / 2^20)}")
