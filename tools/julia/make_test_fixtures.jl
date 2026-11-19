# Generate the fixture files used by the C++ tests (tests/data/), reproducing
# the random inputs that SNaQ.jl's tests create with PhyloCoalSimulations.jl
# and Julia's RNG (simulated gene trees, random networks, starting trees).
#
# usage (from the repository root):
#   julia --project=tools/julia tools/julia/make_test_fixtures.jl
#
# Fixture formats:
#   *.net / *.tre : one Newick network or tree per line
#   *.cf          : first line "nrow", then one row "q1 q2 q3" per quartet
#   *.txt         : documented where written

using SNaQ, PhyloNetworks, PhyloCoalSimulations, Random, StatsBase, Printf, CSV, DataFrames

const SNAQ_DIR = isdir(joinpath(@__DIR__, "..", "..", "SNaQ.jl")) ?
                 joinpath(@__DIR__, "..", "..", "SNaQ.jl") :
                 joinpath(@__DIR__, "..", "..", "..", "SNaQ.jl")
const SNAQTEST = joinpath(SNAQ_DIR, "test")
const OUT = joinpath(@__DIR__, "..", "..", "tests", "data")
mkpath(OUT)
include(joinpath(SNAQTEST, "test_inplace_updates", "misc.jl"))  # generate_tree, generate_net

nw(net) = writenewick(net)
function writelines(file, lines)
    open(joinpath(OUT, file), "w") do f
        for l in lines
            println(f, l)
        end
    end
end
function writecf(file, q::AbstractMatrix)
    open(joinpath(OUT, file), "w") do f
        println(f, size(q, 1))
        for i in 1:size(q, 1)
            println(f, join((repr(q[i, k]) for k in 1:3), " "))
        end
    end
end
qmatrix(q) = [q[j].data[k] for j in eachindex(q), k in 1:3]

# n1.netfile (test_opt_API.jl, test_search_API.jl)
cp(joinpath(SNAQTEST, "test_gradient_opt", "n1.netfile"), joinpath(OUT, "n1.netfile"); force=true)

# --- test_opt_API.jl: "-logPL strictly improves" (with and without ρ) ---
# one line per case: "<rho> <newick>", and q matrices opt_improve_<k>.cf
let lines = String[], ntested = 0, nfail = 0
    while ntested < 100
        Random.seed!(ntested + nfail)
        net = generate_net(7, 1, ntested + nfail)
        while shrink2cycles!(net) || shrink3cycles!(net) continue end
        if net.numhybrids == 0
            nfail += 1
            continue
        end
        gts = simulatecoalescent(net, 100, 1)
        SNaQ.semidirectnetwork!(net)
        q, _ = countquartetsintrees(gts, showprogressbar=false)
        ρ = rand() < 0.25 ? 1.0 : rand()
        ntested += 1
        push!(lines, "$(repr(ρ)) $(nw(net))")
        writecf("opt_improve_$(ntested).cf", qmatrix(q))
    end
    writelines("opt_improve.txt", lines)
end

# --- test_rNNI_inplace.jl ---
let
    Random.seed!(42)
    net = generate_net(15, 2)
    rootonedge!(net, getroot(net).edge[1])
    gts = simulatecoalescent(net, 100, 1)
    q, _ = countquartetsintrees(gts, showprogressbar=false)
    SNaQ.semidirectnetwork!(net)
    writelines("inplace.net", [nw(net)])
    writecf("inplace.cf", qmatrix(q))
end

# --- test_search_API.jl ---
let
    net = readnewick(joinpath(SNAQTEST, "test_gradient_opt", "n1.netfile"))
    tre0 = majortree(net)
    writelines("search_tre0.tre", [nw(tre0)])
    # sets of 3 random starting trees (5 uses) and single starting trees (5 uses)
    writelines("search_t0s.tre", [nw(t) for t in simulatecoalescent(tre0, 15, 1)])
    writelines("search_t0.tre", [nw(t) for t in simulatecoalescent(tre0, 5, 1)])

    # level-2 truth outside the restrictions
    truenet = readnewick("(((t4:0.4,(t8:0.3,t2:0.3):0.1):0.2,(t1:0.4,((t5:0.2,#H2:0.0::0.25):0.2,((t6:0.0,(t3:0.0)#H1:0.0::0.75):0.2)#H2:0.2::0.75):0.0):0.2):0.6,(t7:0.6,#H1:0.0::0.25):0.6);")
    writelines("search_level2_T0.tre", [nw(simulatecoalescent(truenet, 1, 1)[1])])

    # many hybrids: "<n> <h> <seed>" then truth and starting tree
    lines = String[]
    for n in [10, 20], h in [3, 6, 9], seed in 1:10
        truenet = generate_net(n, h, seed)
        T0 = simulatecoalescent(truenet, 1, 1)[1]
        push!(lines, "$n $h $seed")
        push!(lines, nw(truenet))
        push!(lines, nw(T0))
    end
    writelines("search_manyhybrids.txt", lines)

    writelines("gennet_10_2_8.net", [nw(generate_net(10, 2, 8))])
    writelines("gennet_10_2_12.net", [nw(generate_net(10, 2, 12))])

    multifurcation = generate_net(10, 3, 41)
    shrinkedges = sample([E for E in multifurcation.edge if !getchild(E).leaf && !E.hybrid], 5, replace=false)
    for E in shrinkedges
        PhyloNetworks.shrinkedge!(multifurcation, E)
    end
    writelines("multifurcation_net.net", [nw(multifurcation)])
end

# --- test_add_remove_retic.jl: major tree of neth2 ---
let
    neth2 = readnewick("(((((D:1,C:1):1)#H2:1::0.8,((A:1,B:1):1)#H1:1::0.7):1,(#H1:1::0.3,E:1):1):1,(#H2:1::0.2,O:1):1);")
    SNaQ.semidirectnetwork!(neth2)
    writelines("neth2_majortree.tre", [nw(majortree(neth2))])
end

# --- test_propQuartets.jl: lazy DataCF gene trees ---
let
    truetre = readnewick("(((a,b):0.5,(c,d):0.5):0.5,(e,f):0.5);")
    simtre = deepcopy(truetre)
    for E in simtre.edge if getchild(E).leaf E.length = 1.0 end end
    Random.seed!(42)
    writelines("propq_gts.tre", [nw(g) for g in simulatecoalescent(simtre, 1000, 1)])
end

# --- test_lazy_datacf.jl ---
let
    Random.seed!(42)
    tre = readnewick("(t1, t2, t3, t4, t5, t6, t7, t8);")
    for E in tre.edge E.length = 0.0 end
    tre = simulatecoalescent(tre, 1, 1)[1]
    gts = simulatecoalescent(tre, 1000, 1)
    writelines("lazy_tre.tre", [nw(tre)])
    writelines("lazy_gts.tre", [nw(g) for g in gts])
end

# --- test_bugs.jl: BUCKy CF table, as the CF matrix SNaQ uses internally ---
let
    bucky_path = joinpath(SNAQTEST, "..", "examples", "buckyCF.csv")
    cp(bucky_path, joinpath(OUT, "buckyCF.csv"); force=true)
    bucky = readtableCF(bucky_path)
    writelines("bucky_taxa.txt", [join(sort(tiplabels(bucky)), " ")])
    writecf("bucky.cf", SNaQ.gatherCFmatrix(bucky))
end

# --- CF tables generated with PhyloNetworks and SNaQ machinery ---
let
    Random.seed!(1234)
    net = readnewick("(((A:1.0,B:1.0):0.5,(C:1.0,D:1.0):0.5):0.5,(E:1.0,F:1.0):0.5);")
    gts = simulatecoalescent(net, 200, 1)
    q, t = countquartetsintrees(gts; showprogressbar=false)
    df = DataFrame(tablequartetCF(q, t))
    CSV.write(joinpath(OUT, "phylonetworks_tableCF.csv"), df)

    dat = readtableCF(df)
    writecf("phylonetworks_table.cf", SNaQ.gatherCFmatrix(dat))
    writelines("phylonetworks_taxa.txt", [join(sort(tiplabels(dat)), " ")])
end

let
    Random.seed!(5678)
    net = readnewick("((((A:0.8,B:0.8):0.4,C:1.2):0.3,D:1.5):0.2,(E:1.0,F:1.0):0.7);")
    gts = simulatecoalescent(net, 150, 1)
    sub_gts = [deepcopy(gt) for gt in gts]
    for (i, gt) in enumerate(sub_gts)
        if i % 3 == 0
            PhyloNetworks.deleteleaf!(gt, "F")
        end
        if i % 5 == 0
            PhyloNetworks.deleteleaf!(gt, "A")
        end
    end
    q, t = countquartetsintrees(sub_gts; showprogressbar=false)
    df = DataFrame(tablequartetCF(q, t))
    CSV.write(joinpath(OUT, "phylonetworks_subsampled_ngenes.csv"), df)

    dat = readtableCF(df)
    writecf("phylonetworks_subsampled.cf", SNaQ.gatherCFmatrix(dat))
    writelines("phylonetworks_subsampled_taxa.txt", [join(sort(tiplabels(dat)), " ")])
end

println("fixtures written to $(abspath(OUT))")
