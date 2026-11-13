# Generate the benchmark datasets in benchmarks/datasets/ (run once; the
# datasets are committed so that every benchmark run uses identical inputs).
#
# usage (from the repository root):
#   julia --project=tools/julia benchmarks/julia/make_datasets.jl
#
# Each dataset <name>/ contains:
#   truth.net      the true network (from checks.net)
#   genetrees.tre  gene trees simulated from it (PhyloCoalSimulations.jl)
#   start.tre      a starting tree: the major tree of the true network,
#                  with branch lengths of at least 0.1

using SNaQ, PhyloNetworks, PhyloCoalSimulations, Random, CSV

const ROOT = joinpath(@__DIR__, "..", "..")
const OUT = joinpath(ROOT, "benchmarks", "datasets")

# name => (number of taxa, number of hybrids, number of gene trees)
const DATASETS = [
    ("n10_h1", 10, 1, 500),
    ("n15_h1", 15, 1, 500),
    ("n20_h2", 20, 2, 500),
    ("n25_h1", 25, 1, 500),
    ("n25_h3", 25, 3, 500),
    ("n30_h2", 30, 2, 500),
    ("n40_h1", 40, 1, 1000),
    ("n50_h2", 50, 2, 1000),
]

function main()
    lines = readlines(joinpath(ROOT, "checks.net"))
    res = CSV.File(joinpath(ROOT, "validation", "checks_net_comparison.csv"))
    Random.seed!(2026)
    for (name, n, h, ngt) in DATASETS
        cand = [r.line for r in res if r.status == "match" && r.ntaxa == n && r.nhybrids == h]
        isempty(cand) && error("no valid network in checks.net with $n taxa and $h hybrids")
        net = readnewick(lines[cand[1]])
        gts = simulatecoalescent(net, ngt, 1)
        start = majortree(net)
        for e in start.edge
            e.length = max(e.length, 0.1)
        end
        dir = joinpath(OUT, name)
        mkpath(dir)
        write(joinpath(dir, "truth.net"), lines[cand[1]] * "\n")
        open(joinpath(dir, "genetrees.tre"), "w") do f
            for g in gts
                println(f, writenewick(g))
            end
        end
        write(joinpath(dir, "start.tre"), writenewick(start) * "\n")
        println("$name: checks.net line $(cand[1]), $ngt gene trees")
    end
end

main()
