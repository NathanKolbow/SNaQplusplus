# Benchmark numerical parameter optimization: SNaQ.jl's fitnumericalparameters!
# versus CSNaQ (tools/bench_fit.cpp), on the same networks, observed CFs and
# starting values.
#
# usage (from the repository root, after building):
#   julia -t <threads> --project=tools/julia tools/julia/bench_fit.jl build/csnaq_bench_fit cases.txt out.csv [maxeval]
# where each line of cases.txt is "<network newick file> <gene trees file>".

using SNaQ, PhyloNetworks, Printf

function juliafit(netfile, gtfile, maxeval; repeats=3)
    net0 = readnewick(readlines(netfile)[1])
    trees = readmultinewick(gtfile)
    d = redirect_stdout(devnull) do
        DataCF(trees)
    end
    q = SNaQ.gatherCFmatrix(d)
    best_fit = Inf; best_eq = Inf; score = NaN; np = 0; net = nothing
    for r in 1:repeats
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
        best_fit = min(best_fit, tfit); best_eq = min(best_eq, teq)
        np = SNaQ.countoptimizationparams(net)
    end
    return (; ntaxa=net.numtaxa, nhyb=net.numhybrids, nq=size(q, 1), np, score, best_fit, best_eq)
end

function main(args)
    bench, cases, outcsv = args[1], args[2], args[3]
    maxeval = length(args) >= 4 ? parse(Int, args[4]) : 1000
    nthreads = Threads.nthreads()
    open(outcsv, "w") do out
        println(out, "network,ntaxa,nhybrids,nquartets,nparams,maxeval,threads,julia_score,cpp_score,abs_score_diff,julia_fit_s,cpp_fit_s,speedup_fit,julia_eqns_s,cpp_eqns_s,speedup_eqns")
        for line in readlines(cases)
            isempty(strip(line)) && continue
            netfile, gtfile = split(line)
            j = juliafit(netfile, gtfile, maxeval)
            c = split(readchomp(`$bench $netfile $gtfile $maxeval $nthreads`))
            cs = parse(Float64, c[5]); cfit = parse(Float64, c[6]); ceq = parse(Float64, c[7])
            row = @sprintf("%s,%d,%d,%d,%d,%d,%d,%.12f,%.12f,%.3e,%.4f,%.4f,%.1f,%.4f,%.4f,%.1f",
                basename(netfile), j.ntaxa, j.nhyb, j.nq, j.np, maxeval, nthreads, j.score, cs, abs(j.score - cs),
                j.best_fit, cfit, j.best_fit / cfit, j.best_eq, ceq, j.best_eq / ceq)
            println(out, row); flush(out)
            println(row)
        end
    end
end

main(ARGS)
