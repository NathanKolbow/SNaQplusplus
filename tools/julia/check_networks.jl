# Compare SNaQ++ against SNaQ.jl on every network of a Newick file:
# expected quartet CFs (ρ = 0 and ρ = 0.5), SNaQ scores and gradients
# (ρ = 0 with and without the tree-quartet batch, and ρ = 0.5).
#
# usage (from the repository root, after building):
#   julia -t 8 --project=tools/julia tools/julia/check_networks.jl checks.net build/snaqpp_check out.csv [cap] [first] [last]
#
# For networks with more than `cap` quartets (default 20000), a fixed,
# evenly spaced subset of `cap` quartets is compared (the same subset on
# both sides): SNaQ.jl takes hours per network on all C(100,4) quartets.
# Synthetic observed CFs are a deterministic function of the quartet rank,
# defined identically in tools/snaqpp_check.cpp.

using SNaQ, PhyloNetworks, Printf

function synthetic_obs(r::Int)  # r: 0-based rank
    a = 1 + ((r * 7) % 11); b = 1 + ((r * 13 + 5) % 7); c = 1 + ((r * 3 + 1) % 5)
    r % 4 == 0 && (c = 0)
    r % 9 == 0 && (a = 0)
    s = a + b + c
    return (a / s, b / s, c / s)
end

function quartet_subset(nq::Int, cap::Int)  # 1-based ranks
    (cap <= 0 || nq <= cap) && return collect(1:nq)
    return [div(i * nq, cap) + 1 for i in 0:(cap-1)]
end

function readcpp(path)
    open(path) do io
        nq = read(io, Int64); np = read(io, Int64)
        s0 = read(io, Float64); sb = read(io, Float64); s05 = read(io, Float64)
        rd(n) = [read(io, Float64) for _ in 1:n]
        g0 = rd(np); gb = rd(np); g05 = rd(np)
        e0 = reshape(rd(3nq), 3, nq)'; e05 = reshape(rd(3nq), 3, nq)'
        return (; nq, np, s0, sb, s05, g0, gb, g05, e0, e05)
    end
end

function juliavalues(line::String, cap::Int)
    net = readnewick(line)
    N = SNaQ.deepcopynetwork(net)
    SNaQ.semidirectnetwork!(N)
    nqall = SNaQ.nchoose4taxalength(N)
    ranks = quartet_subset(nqall, cap)
    eqns, _, params, _, _ = SNaQ.findquartetequations(N, ranks)
    nq = length(ranks)
    obs = zeros(nq, 3)
    for (j, r) in enumerate(ranks)
        obs[j, :] .= synthetic_obs(r - 1)
    end
    np = length(params)
    g0 = zeros(np); s0 = SNaQ.computelossandgradient!(eqns, params, g0, obs, Inf)
    g05 = zeros(np); s05 = SNaQ.computelossandgradient!(eqns, params, g05, obs, SNaQ.rhotoalpha(0.5))
    batch = SNaQ.treequartetbatch(eqns, obs)
    gb = zeros(np); sb = SNaQ.computelossandgradient!(batch, eqns, params, gb, obs, Inf)
    ecf(ρ) = begin
        M = zeros(nq, 3)
        for j in 1:nq
            e1, e2 = SNaQ.computeexpectedCF(eqns[j], params, ρ)
            M[j, :] .= (e1, e2, 1 - e1 - e2)
        end
        M
    end
    return (; ntaxa=N.numtaxa, nhyb=N.numhybrids, nq, np, s0, sb, s05, g0, gb, g05, e0=ecf(0.0), e05=ecf(0.5))
end

maxdiff(a, b) = isempty(a) ? 0.0 : maximum(abs.(a .- b))

function main(args)
    netfile, checker, outcsv = args[1], args[2], args[3]
    cap = length(args) >= 4 ? parse(Int, args[4]) : 20000
    lines = readlines(netfile)
    first = length(args) >= 5 ? parse(Int, args[5]) : 1
    last = length(args) >= 6 ? parse(Int, args[6]) : length(lines)
    tmp = tempname()
    nthreads = Threads.nthreads()
    worst = 0.0
    nbad = 0
    open(outcsv, "w") do out
        println(out, "line,ntaxa,nhybrids,nquartets,nparams,status,max_ecf_rho0,max_ecf_rho05,diff_score_rho0,diff_score_batch,diff_score_rho05,max_grad_rho0,max_grad_batch,max_grad_rho05")
        for i in first:last
            jerr = nothing; jv = nothing
            try
                jv = juliavalues(lines[i], cap)
            catch e
                jerr = sprint(showerror, e)
            end
            ok = success(pipeline(`$checker $netfile $i $tmp $cap $nthreads`, stderr=devnull))
            status = ""
            row = ""
            if jerr !== nothing || !ok
                status = (jerr !== nothing && !ok) ? "both_error" :
                         (jerr !== nothing ? "julia_error_only" : "cpp_error_only")
                status == "both_error" || (nbad += 1)
                row = @sprintf("%d,,,,,%s,,,,,,,,", i, status)
            else
                cv = readcpp(tmp)
                if cv.nq != jv.nq || cv.np != jv.np
                    status = "size_mismatch"; nbad += 1
                    row = @sprintf("%d,%d,%d,%d,%d,%s,,,,,,,,", i, jv.ntaxa, jv.nhyb, jv.nq, jv.np, status)
                else
                    d = (maxdiff(jv.e0, cv.e0), maxdiff(jv.e05, cv.e05), abs(jv.s0 - cv.s0), abs(jv.sb - cv.sb),
                         abs(jv.s05 - cv.s05), maxdiff(jv.g0, cv.g0), maxdiff(jv.gb, cv.gb), maxdiff(jv.g05, cv.g05))
                    m = maximum(d)
                    worst = max(worst, m)
                    status = m <= 1e-8 ? "match" : "MISMATCH"
                    status == "match" || (nbad += 1)
                    row = @sprintf("%d,%d,%d,%d,%d,%s,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e", i, jv.ntaxa,
                                   jv.nhyb, jv.nq, jv.np, status, d...)
                end
            end
            println(out, row); flush(out)
            i % 25 == 0 && println("checked network $i ($(nbad) problem(s) so far, worst diff $(worst))")
        end
    end
    isfile(tmp) && rm(tmp)
    println("done: $(nbad) network(s) with a mismatch; worst difference among matched networks: $(worst)")
    return nbad
end

exit(main(ARGS) == 0 ? 0 : 1)
