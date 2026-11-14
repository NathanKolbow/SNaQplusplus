# Test of the restriction flags of the snaq program (--max_level, --galled_tree,
# --galled_network, --rooted_tree_child, --weakly_tree_child,
# --strongly_tree_child): runs snaq with each of them, reads the networks it
# wrote (the best network and the network of every run, from <filename>.out
# and <filename>.networks) and checks every one of them against the same
# restriction with SNaQ.jl's own functions (SNaQ.restrictionset).
#
# Each network is checked as SNaQ.jl's search checks its candidates: on a copy
# made with SNaQ.deepcopynetwork, after SNaQ.semidirectnetwork!.
#
# A run without restrictions shows how many of its networks each restriction
# would reject: a restriction that rejects none of them is not exercised by
# these data (reported, not an error).
#
# usage (from the repository root, after building):
#   julia --project=tools/julia tools/julia/check_restrictions.jl <snaq binary> [dataset directory] [hmax] [runs] [Nfail]
# The dataset directory holds genetrees.tre and start.tre (default:
# benchmarks/datasets/n20_h2). Exits with 1 if a network violates the
# restriction it was searched under, or if a run fails.

using SNaQ, PhyloNetworks, Printf

const binary = abspath(ARGS[1])
const dataset = length(ARGS) > 1 ? ARGS[2] : joinpath(@__DIR__, "..", "..", "benchmarks", "datasets", "n20_h2")
const hmax = length(ARGS) > 2 ? parse(Int, ARGS[3]) : 4
const runs = length(ARGS) > 3 ? parse(Int, ARGS[4]) : 16
const nfail = length(ARGS) > 4 ? parse(Int, ARGS[5]) : 25

# (name, snaq flags, keyword arguments of SNaQ.restrictionset)
const cases = [
    ("max_level 0", ["--max_level", "0"], (max_level = 0,)),
    ("max_level 1", ["--max_level", "1"], (max_level = 1,)),
    ("max_level 2", ["--max_level", "2"], (max_level = 2,)),
    ("galled_tree", ["--galled_tree"], (galled_tree = true,)),
    ("galled_network", ["--galled_network"], (galled_network = true,)),
    ("rooted_tree_child", ["--rooted_tree_child"], (rooted_tree_child = true,)),
    ("weakly_tree_child", ["--weakly_tree_child"], (weakly_tree_child = true,)),
    ("strongly_tree_child", ["--strongly_tree_child"], (strongly_tree_child = true,)),
    ("max_level 2 + galled_network", ["--max_level", "2", "--galled_network"], (max_level = 2, galled_network = true)),
    ("max_level 2 + strongly_tree_child", ["--max_level=2", "--strongly_tree_child=true"],
     (max_level = 2, strongly_tree_child = true)),
    ("galled_network + weakly_tree_child", ["--galled_network", "--weakly_tree_child"],
     (galled_network = true, weakly_tree_child = true)),
    # galled_tree takes precedence over galled_network, rooted over weakly over strongly
    ("galled_tree + galled_network", ["--galled_tree", "--galled_network"], (galled_tree = true, galled_network = true)),
    ("rooted + strongly tree-child", ["--rooted_tree_child", "--strongly_tree_child"],
     (rooted_tree_child = true, strongly_tree_child = true)),
    # hyphenated spellings of the flags
    ("max-level 1 (hyphenated)", ["--max-level", "1"], (max_level = 1,)),
    ("galled-network (hyphenated)", ["--galled-network"], (galled_network = true,)),
    ("strongly-tree-child (hyphenated)", ["--strongly-tree-child"], (strongly_tree_child = true,)),
]

# The check prescribed for SNaQ.jl: copy, semidirect, then test the copy.
function meets(net::HybridNetwork, restriction)::Bool
    n = SNaQ.deepcopynetwork(net)
    SNaQ.semidirectnetwork!(n)
    return restriction(n)
end

# The networks of a snaq output file: every line that starts with a Newick string.
function networksof(file::String)
    nets = HybridNetwork[]
    for line in eachline(file)
        m = match(r"^\s*(\(.*?;)", line)
        m === nothing && continue
        push!(nets, readnewick(String(m.captures[1])))
    end
    return nets
end

function runsnaq(flags::Vector{String}, outroot::String, seed::Int)
    cmd = `$binary -i $(joinpath(dataset, "genetrees.tre")) --start-tree $(joinpath(dataset, "start.tre")) -o $outroot
           --hmax $hmax --runs $runs --Nfail $nfail --seed $seed --propQuartets 0.2 --propQuartetsFinal 0.2 $flags`
    err = IOBuffer()
    ok = success(pipeline(cmd; stdout = devnull, stderr = err))
    ok || println("  snaq failed: ", String(take!(err)))
    return ok
end

function main()
    tmp = mktempdir()
    failures = 0

    # control: no restriction
    println("snaq: ", binary, "\ndata: ", dataset, " (hmax ", hmax, ", ", runs, " runs, Nfail ", nfail, ")\n")
    control = HybridNetwork[]
    if runsnaq(String[], joinpath(tmp, "control"), 1)
        control = vcat(networksof(joinpath(tmp, "control.out")), networksof(joinpath(tmp, "control.networks")))
    else
        failures += 1
    end
    @printf("%-36s %9s %9s %10s   %s\n", "restriction", "networks", "violate", "max hybrids", "rejected among the networks found without restriction")
    for (name, flags, kw) in cases
        restriction = SNaQ.restrictionset(; kw...)
        outroot = joinpath(tmp, replace(name, r"[^A-Za-z0-9]+" => "_"))
        if !runsnaq(flags, outroot, 1)
            failures += 1
            continue
        end
        # the best network (first line of .out), the runs listed in .out, and .networks
        nets = vcat(networksof(outroot * ".out"), networksof(outroot * ".networks"))
        if isempty(nets)
            println(name, ": no network in the output files")
            failures += 1
            continue
        end
        bad = count(n -> !meets(n, restriction), nets)
        rejected = count(n -> !meets(n, restriction), control)
        @printf("%-36s %9d %9d %10d   %d of %d%s\n", name, length(nets), bad, maximum(n.numhybrids for n in nets),
                rejected, length(control), rejected == 0 ? "  (not exercised by these data)" : "")
        failures += bad
    end
    rm(tmp; recursive = true, force = true)
    if failures > 0
        println("\nFAILED: ", failures, " violation(s) or failed run(s)")
        exit(1)
    end
    println("\nall networks meet the restriction they were searched under")
end

main()
