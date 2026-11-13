# One end-to-end benchmark case with SNaQ.jl's snaq! (called by benchmarks/bench.py).
#
# usage: julia --project=tools/julia run_snaq_case.jl <genetrees> <start tree> <hmax> <runs> <Nfail>
#            <seed> <propQuartets> <propQuartetsFinal> <lazy: true|false> <output root>
#
# Prints a line "RESULT {json}" with
#   snaq_s     wall time of the snaq! call (after a warm-up call, so Julia's
#              compilation is not included)
#   obscf_s    wall time of building the DataCF from the gene trees
#              (countquartetsintrees, PhyloNetworks' fast algorithm, for a
#              non-lazy DataCF; DataCF(trees; lazy=true) otherwise)
#   peak_rss_mb / base_rss_mb  peak resident set size of the process, and
#              after loading the packages
using SNaQ, PhyloNetworks, DataFrames

base = Sys.maxrss()
gtfile, startfile = ARGS[1], ARGS[2]
hmax, runs, nfail, seed = parse.(Int, ARGS[3:6])
pq, pqf = parse.(Float64, ARGS[7:8])
lazy = ARGS[9] == "true"
outroot = ARGS[10]

function makedata(trees, lazy)
    lazy && return DataCF(trees; lazy=true)
    q, t = countquartetsintrees(trees; showprogressbar=false)
    return readtableCF(DataFrame(tablequartetCF(q, t)))
end

trees = readmultinewick(gtfile)
start = readnewick(readlines(startfile)[1])

# warm-up (compilation) on a small problem, with the same kind of DataCF
wt = readnewick("((a:1,b:1):1,(c:1,d:1):1,(e:1,f:1):1);")
wgts = [readnewick("((a,b),(c,d),(e,f));"), readnewick("((a,c),(b,d),(e,f));"), readnewick("((a,b),(c,e),(d,f));")]
redirect_stdout(devnull) do
    wd = makedata(wgts, lazy)
    snaq!(wt, wd; hmax=1, runs=1, Nfail=2, filename="", seed=1,
          propQuartets=lazy ? 0.9 : 1.0, propQuartetsFinal=lazy ? 0.9 : 1.0)
end

tcf = @elapsed d = redirect_stdout(devnull) do
    makedata(trees, lazy)
end
t = @elapsed net = snaq!(start, d; hmax=hmax, runs=runs, Nfail=nfail, seed=seed,
                         propQuartets=pq, propQuartetsFinal=pqf, filename=outroot)
println("RESULT {\"snaq_s\": $t, \"obscf_s\": $tcf, \"reported_score\": $(SNaQscore(net)), ",
        "\"peak_rss_mb\": $(Sys.maxrss() / 2^20), \"base_rss_mb\": $(base / 2^20), ",
        "\"julia_threads\": $(Threads.nthreads())}")
