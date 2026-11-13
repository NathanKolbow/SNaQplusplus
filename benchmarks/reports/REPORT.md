# CSNaQ benchmark report

Generated 2026-11-13T16:42:36-06:00 by `benchmarks/bench.py report` from `benchmarks/history/results.jsonl`. Baseline: **SNaQ.jl**.

Time is the wall time of the whole `snaq` process for C++. For SNaQ.jl it is the time to build the DataCF plus the `snaq!` call, after a warm-up call, so Julia's start-up and compilation are not included. Peak memory is the peak resident set size of the process; for SNaQ.jl it includes the Julia runtime and loaded packages (see `base_rss_mb` in the history). HWCD is the hardwired cluster distance between the best network and the true network (0: same topology). The score on all quartets is the SNaQ score of the best network recomputed by `csnaq_eval` on the observed CFs of all quartets, which makes runs that used `propQuartets` < 1 comparable. "Refitted to all quartets" is the score of the best network's topology after its parameters are fitted to all quartets, as is done for the true network: it compares topologies, whatever quartets the reported parameters were estimated from. After `bench.py rescore`, HWCD and the scores on all quartets of every version come from the same `csnaq_eval`.

## Suite `standard`

End-to-end snaq runs on simulated data from 10 to 50 taxa and 1 to 3 reticulations, 100 runs each. Each case uses about 1000 quartets: propQuartets = propQuartetsFinal = min(1, 1000 / binomial(ntaxa, 4)). Cases with 20 or more taxa use --lazy (observed CFs computed on demand).

Latest run of each version: **SNaQ.jl**: run `20261113T141603_SNaQ.jl`, **exp-C-newton**: run `20261113T153817_exp-C-newton`, **exp-C-newton-hybrids**: run `20261113T154448_exp-C-newton-hybrids`, **exp-D-eqbuilder**: run `20261113T163336_exp-D-eqbuilder`.

![Runtime](standard_time.svg)

![Peak memory](standard_memory.svg)

![Distance of the best network to the truth](standard_hwcd.svg)

![Score of the best network on all quartets, minus the true network's](standard_score.svg)

![Score of the best topology refitted to all quartets, minus the true network's](standard_refit.svg)

![Speedup](standard_speedup.svg)

| case | version | time (s) | peak MB | best score | median run score | score on all quartets | refitted to all quartets | truth score on all quartets | HWCD | hybrids |
|---|---|---|---|---|---|---|---|---|---|---|
| n10_h1 | SNaQ.jl | 24.9132 | 1143.64 | -0.223569 | -0.270387 | -0.233167 | -0.226623 | -0.22524 | 4 | 1 |
| n10_h1 | exp-C-newton | 0.933821 | 8.53125 | -0.22524 | -0.22673 | -0.22524 | -0.22524 | -0.22524 | 0 | 1 |
| n10_h1 | exp-C-newton-hybrids | 0.739359 | 8.73438 | -0.22524 | -0.22673 | -0.22524 | -0.22524 | -0.22524 | 0 | 1 |
| n10_h1 | exp-D-eqbuilder | 0.791147 | 8.73438 | -0.22524 | -0.22673 | -0.22524 | -0.22524 | -0.22524 | 0 | 1 |
| n15_h1 | SNaQ.jl | 208.937 | 1090.34 | -1.43115 | -3.4746 | -1.98535 | -1.9832 | -2.20124 | 6 | 1 |
| n15_h1 | exp-C-newton | 10.1716 | 9.46875 | -1.40239 | -3.09541 | -1.98648 | -1.9832 | -2.20124 | 6 | 1 |
| n15_h1 | exp-C-newton-hybrids | 10.5675 | 9.17188 | -1.40252 | -2.68042 | -1.98414 | -1.9832 | -2.20124 | 6 | 1 |
| n15_h1 | exp-D-eqbuilder | 4.80601 | 9.46875 | -1.40241 | -2.68042 | -1.9861 | -1.9832 | -2.20124 | 6 | 1 |
| n20_h2 | SNaQ.jl | 180.818 | 1081.59 | -1.21354 | -7.81406 | -6.37721 | -6.32935 | -6.55637 | 12 | 2 |
| n20_h2 | exp-C-newton | 22.0962 | 13.9531 | -1.24589 | -1.40744 | -6.4881 | -6.46009 | -6.55637 | 10 | 2 |
| n20_h2 | exp-C-newton-hybrids | 23.6029 | 13.9844 | -1.21689 | -1.3512 | -6.36573 | -6.33145 | -6.55637 | 10 | 2 |
| n20_h2 | exp-D-eqbuilder | 9.93971 | 14.4531 | -1.21689 | -1.3512 | -6.36573 | -6.33145 | -6.55637 | 10 | 2 |
| n25_h1 | SNaQ.jl | 186.532 | 1184.95 | -1.38173 | -1.56734 | -19.4624 | -19.2281 | -20.1041 | 20 | 1 |
| n25_h1 | exp-C-newton | 31.5951 | 13.8281 | -1.34231 | -1.5081 | -18.2185 | -18.088 | -20.1041 | 20 | 1 |
| n25_h1 | exp-C-newton-hybrids | 30.5803 | 13.625 | -1.35522 | -1.5034 | -18.1914 | -18.1192 | -20.1041 | 22 | 1 |
| n25_h1 | exp-D-eqbuilder | 14.9908 | 14.1719 | -1.35522 | -1.5034 | -18.1914 | -18.1192 | -20.1041 | 22 | 1 |
| n25_h3 | SNaQ.jl | 1166.13 | 968.203 | -1.57573 | -4.62509 | -21.6193 | -21.2298 | -20.1142 | 20 | 3 |
| n25_h3 | exp-C-newton | 189.585 | 14.4062 | -1.48677 | -3.45363 | -20.0842 | -19.8409 | -20.1142 | 10 | 3 |
| n25_h3 | exp-C-newton-hybrids | 192.423 | 15.6719 | -1.56264 | -3.26133 | -20.753 | -20.5659 | -20.1142 | 12 | 3 |
| n25_h3 | exp-D-eqbuilder | 85.3607 | 15.25 | -1.56264 | -3.32811 | -20.753 | -20.5659 | -20.1142 | 12 | 3 |
| n30_h2 | SNaQ.jl | 312.185 | 966.828 | -1.64993 | -2.28227 | -49.1109 | -48.5461 | -48.9468 | 8 | 2 |
| n30_h2 | exp-C-newton | 52.3083 | 14.7656 | -1.5993 | -1.86896 | -46.358 | -45.7432 | -48.9468 | 16 | 2 |
| n30_h2 | exp-C-newton-hybrids | 50.2599 | 14.9375 | -1.56359 | -1.8014 | -47.1689 | -46.7394 | -48.9468 | 12 | 2 |
| n30_h2 | exp-D-eqbuilder | 22.2814 | 14.875 | -1.59274 | -1.83441 | -49.6484 | -48.9063 | -48.9468 | 12 | 2 |
| n40_h1 | SNaQ.jl | 210.55 | 1209.11 | -0.826667 | -1.09699 | -83.5046 | -82.0908 | -81.2316 | 2 | 1 |
| n40_h1 | exp-C-newton | 24.2371 | 34.4531 | -0.821061 | -1.07701 | -81.7989 | -80.2382 | -81.2316 | 22 | 1 |
| n40_h1 | exp-C-newton-hybrids | 26.4365 | 33.7969 | -0.823797 | -1.08055 | -81.5232 | -80.4679 | -81.2316 | 18 | 1 |
| n40_h1 | exp-D-eqbuilder | 15.4492 | 33.8281 | -0.823797 | -1.08055 | -81.5232 | -80.4679 | -81.2316 | 18 | 1 |
| n50_h2 | SNaQ.jl | 389.614 | 1131.92 | -1.15819 | -18.212 | -266.663 | -256.563 | -197.667 | 14 | 2 |
| n50_h2 | exp-C-newton | 56.9988 | 36.5 | -0.830043 | -8.44807 | -203.798 | -200.408 | -197.667 | 10 | 2 |
| n50_h2 | exp-C-newton-hybrids | 59.1139 | 36.5156 | -0.844858 | -7.23005 | -207.139 | -202.111 | -197.667 | 16 | 2 |
| n50_h2 | exp-D-eqbuilder | 29.6688 | 36.25 | -0.904829 | -7.20874 | -215.027 | -209.914 | -197.667 | 10 | 2 |

![History](standard_history.svg)

## Suite `fit`

Numerical parameter optimization of a fixed network (the true network) on the observed CFs of all quartets: fitnumericalparameters! with maxeval = 1000.

Latest run of each version: **SNaQ.jl**: run `20261113T150153_SNaQ.jl`, **exp-A-classes**: run `20261113T151955_exp-A-classes`, **exp-C-newton-hybrids**: run `20261113T155124_exp-C-newton-hybrids`, **exp-D-eqbuilder**: run `20261113T163642_exp-D-eqbuilder`.

![Parameter optimization time](fit_fit.svg)

![Time to build the quartet equations](fit_equations.svg)

![Peak memory](fit_memory.svg)

![Speedup](fit_speedup.svg)

| case | version | fit (s) | equations (s) | peak MB | score | quartets | parameters |
|---|---|---|---|---|---|---|---|
| fit_n10_h1 | SNaQ.jl | 0.00596933 | 0.00514454 | 692.875 | -0.22524 | 210 | 11 |
| fit_n10_h1 | exp-A-classes | 0.000572125 | 0.000295208 | 6.23438 | -0.22524 | 210 | 11 |
| fit_n10_h1 | exp-C-newton-hybrids | 0.000586875 | 0.000292666 | 6.29688 | -0.22524 | 210 | 11 |
| fit_n10_h1 | exp-D-eqbuilder | 0.000624542 | 5.0833e-05 | 6.04688 | -0.22524 | 210 | 11 |
| fit_n15_h1 | SNaQ.jl | 0.0775856 | 0.0265386 | 715.984 | -2.20124 | 1365 | 16 |
| fit_n15_h1 | exp-A-classes | 0.000678875 | 0.0007845 | 6.6875 | -2.20124 | 1365 | 16 |
| fit_n15_h1 | exp-C-newton-hybrids | 0.0004865 | 0.000832958 | 6.6875 | -2.20124 | 1365 | 16 |
| fit_n15_h1 | exp-D-eqbuilder | 0.000451666 | 0.000164459 | 6.5 | -2.20124 | 1365 | 16 |
| fit_n20_h2 | SNaQ.jl | 0.290687 | 0.0717673 | 749.812 | -6.55637 | 4845 | 23 |
| fit_n20_h2 | exp-A-classes | 0.00387 | 0.00274704 | 10.2344 | -6.55637 | 4845 | 23 |
| fit_n20_h2 | exp-C-newton-hybrids | 0.00105467 | 0.00250167 | 10.3125 | -6.55637 | 4845 | 23 |
| fit_n20_h2 | exp-D-eqbuilder | 0.00111958 | 0.000610583 | 10.1875 | -6.55637 | 4845 | 23 |
| fit_n25_h1 | SNaQ.jl | 0.776733 | 0.154207 | 711.922 | -20.1041 | 12650 | 24 |
| fit_n25_h1 | exp-A-classes | 0.00852912 | 0.00354917 | 10.9219 | -20.1041 | 12650 | 24 |
| fit_n25_h1 | exp-C-newton-hybrids | 0.00140842 | 0.00363904 | 10.9375 | -20.1041 | 12650 | 24 |
| fit_n25_h1 | exp-D-eqbuilder | 0.00131363 | 0.00107692 | 10.8281 | -20.1041 | 12650 | 24 |
| fit_n25_h3 | SNaQ.jl | 4.43912 | 1.21888 | 859.797 | -27.5921 | 12650 | 32 |
| fit_n25_h3 | exp-A-classes | 0.0237497 | 0.0155119 | 11.5625 | -20.1142 | 12650 | 32 |
| fit_n25_h3 | exp-C-newton-hybrids | 0.00838388 | 0.016571 | 11.6719 | -20.1142 | 12650 | 32 |
| fit_n25_h3 | exp-D-eqbuilder | 0.00669133 | 0.00156421 | 11.3281 | -20.1142 | 12650 | 32 |
| fit_n30_h2 | SNaQ.jl | 1.05106 | 0.51477 | 924.594 | -48.9468 | 27405 | 31 |
| fit_n30_h2 | exp-A-classes | 0.00426458 | 0.00480363 | 11.6406 | -48.9468 | 27405 | 31 |
| fit_n30_h2 | exp-C-newton-hybrids | 0.00257471 | 0.00436304 | 11.6406 | -48.9468 | 27405 | 31 |
| fit_n30_h2 | exp-D-eqbuilder | 0.00209942 | 0.00128379 | 11.4375 | -48.9468 | 27405 | 31 |
| fit_n40_h1 | SNaQ.jl | 1.918 | 1.18424 | 1056.75 | -81.2316 | 91390 | 39 |
| fit_n40_h1 | exp-A-classes | 0.0112363 | 0.00799783 | 35.4688 | -81.2316 | 91390 | 39 |
| fit_n40_h1 | exp-C-newton-hybrids | 0.0141902 | 0.00759613 | 35.6719 | -81.2316 | 91390 | 39 |
| fit_n40_h1 | exp-D-eqbuilder | 0.0117992 | 0.00340462 | 35.4531 | -81.2316 | 91390 | 39 |
| fit_n50_h2 | SNaQ.jl | 27.366 | 8.6391 | 1184.75 | -197.667 | 230300 | 53 |
| fit_n50_h2 | exp-A-classes | 0.065732 | 0.0309066 | 43.875 | -197.667 | 230300 | 53 |
| fit_n50_h2 | exp-C-newton-hybrids | 0.0242798 | 0.0320116 | 44 | -197.667 | 230300 | 53 |
| fit_n50_h2 | exp-D-eqbuilder | 0.0259122 | 0.00968175 | 43.6719 | -197.667 | 230300 | 53 |

![History](fit_history.svg)

