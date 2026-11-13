# End-to-end benchmark of snaq! (SNaQ.jl) versus the snaq program (CSNaQ):
# wall time, SNaQ score and peak memory, each case in a fresh process.
#
# usage (from the repository root, after building):
#   julia --project=tools/julia tools/julia/bench_snaq.jl <snaq binary> <cases file> <out.csv> [threads]
#
# Each line of the cases file is:
#   <name> <gene trees file> <start tree file> <hmax> <runs> <Nfail>
#
# Julia: `snaq!(start, DataCF(genetrees); ...)` is timed after a small warm-up
# call (so compilation is not counted); time to compute the observed CFs is
# reported separately. Peak memory is the peak resident set size of the
# process (Sys.maxrss()); the resident set size after loading the packages
# is reported as `julia_base_MB`.
# C++: the whole `snaq` process is timed (it includes reading the gene trees
# and computing the observed CFs); its peak resident set size is measured with
# /usr/bin/time.

using Printf

const JULIA_CASE = raw"""
using SNaQ, PhyloNetworks
base = Sys.maxrss()
gtfile, startfile, hmax, runs, nfail, outroot = ARGS[1], ARGS[2], parse(Int, ARGS[3]), parse(Int, ARGS[4]), parse(Int, ARGS[5]), ARGS[6]
trees = readmultinewick(gtfile)
start = readnewick(readlines(startfile)[1])
tcf = @elapsed d = redirect_stdout(devnull) do; DataCF(trees); end
# warm-up on a 6-taxon problem: compiles snaq! and the functions it calls
wt = readnewick("((a:1,b:1):1,(c:1,d:1):1,(e:1,f:1):1);")
wgts = [readnewick("((a,b),(c,d),(e,f));"), readnewick("((a,c),(b,d),(e,f));"), readnewick("((a,b),(c,e),(d,f));")]
wd = redirect_stdout(devnull) do; DataCF(wgts); end
redirect_stdout(devnull) do; snaq!(wt, wd; hmax=1, runs=1, Nfail=2, filename="", seed=1); end
t = @elapsed net = snaq!(start, d; hmax=hmax, runs=runs, Nfail=nfail, seed=42, filename=outroot)
println("RESULT ", t, " ", tcf, " ", SNaQscore(net), " ", Sys.maxrss() / 2^20, " ", base / 2^20, " ", Threads.nthreads())
"""

function cpp_run(binary, gt, start, hmax, runs, nfail, outroot, threads)
    timecmd = Sys.isapple() ? `/usr/bin/time -l` : `/usr/bin/time -v`
    cmd = `$timecmd $binary -i $gt --start-tree $start -o $outroot --hmax $hmax --runs $runs --Nfail $nfail --seed 42 --threads $threads`
    out = IOBuffer(); err = IOBuffer()
    t = @elapsed run(pipeline(cmd, stdout=out, stderr=err))
    o = String(take!(out)); e = String(take!(err))
    score = parse(Float64, match(r"SNaQ score = (\S+)", o).captures[1])
    rss = if Sys.apple()
        parse(Float64, match(r"(\d+)\s+maximum resident set size", e).captures[1]) / 2^20
    else
        parse(Float64, match(r"Maximum resident set size \(kbytes\): (\d+)", e).captures[1]) / 2^10
    end
    return t, score, rss
end

function main(args)
    binary, casesfile, outcsv = args[1], args[2], args[3]
    threads = length(args) >= 4 ? parse(Int, args[4]) : 1
    scratch = mktempdir()
    juliascript = joinpath(scratch, "case.jl")
    write(juliascript, JULIA_CASE)
    project = dirname(@__FILE__)
    open(outcsv, "w") do out
        println(out, "case,hmax,runs,Nfail,threads,julia_snaq_s,julia_obsCF_s,cpp_total_s,speedup,julia_score,cpp_score,julia_peak_MB,julia_base_MB,cpp_peak_MB")
        for line in readlines(casesfile)
            isempty(strip(line)) && continue
            name, gt, start, hmax, runs, nfail = split(line)
            jout = read(`$(Base.julia_cmd()) -t $threads --project=$project $juliascript $gt $start $hmax $runs $nfail $(joinpath(scratch, "jl_" * name))`, String)
            r = split(match(r"RESULT (.*)", jout).captures[1])
            jt, jcf, js, jrss, jbase = parse.(Float64, r[1:5])
            ct, cs, crss = cpp_run(binary, gt, start, hmax, runs, nfail, joinpath(scratch, "cpp_" * name), threads)
            row = @sprintf("%s,%s,%s,%s,%d,%.2f,%.2f,%.2f,%.1f,%.10f,%.10f,%.1f,%.1f,%.1f", name, hmax, runs, nfail,
                           threads, jt, jcf, ct, jt / ct, js, cs, jrss, jbase, crss)
            println(out, row); flush(out)
            println(row)
        end
    end
end

main(ARGS)
