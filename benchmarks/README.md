# Benchmarks

A fixed set of benchmark suites, a running history of their results, and figures that compare
versions. Use it to measure every experimental change against SNaQ.jl and against earlier
versions of CSNaQ: runtime, peak memory and accuracy.

**Requirements:**
- Python 3.8 or later, standard library only.
- The C++ tools, built with `make tools`, which also runs the `-DCSNAQ_BUILD_TOOLS=ON` build.
- Julia with the environment in `tools/julia`, only for benchmarking SNaQ.jl. See
  `tools/julia/setup.jl`.

## Quick start

```sh
make bench-quick LABEL="my-idea"              # a few minutes: does the change help?
make bench LABEL="my-idea"                    # the standard suite
make bench-fit LABEL="my-idea"                # parameter optimization only
make bench-report                             # figures + benchmarks/reports/REPORT.md
```

To benchmark experimental `snaq` flags, pass them with `ARGS`. Use a different `LABEL` for
each variant:

```sh
make bench LABEL="new-moves" ARGS="--my-experimental-flag 3"
```

`THREADS=8` sets the number of threads; the default is 1. Results are appended to the
history, so the report compares each version with the baseline and draws the trend over
all recorded runs.

The same commands with the script directly:

```sh
python3 benchmarks/bench.py run --suite standard --label my-idea [--binary path/to/snaq] [--extra-args "..."] [--threads N]
python3 benchmarks/bench.py run --suite fit --label my-idea [--fit-binary build/csnaq_bench_fit]
python3 benchmarks/bench.py run-julia --suite standard           # SNaQ.jl, recorded with label "SNaQ.jl"
python3 benchmarks/bench.py report [--baseline SNaQ.jl] [--labels SNaQ.jl,cpp-v1,my-idea]
python3 benchmarks/bench.py list
python3 benchmarks/bench.py rescore                              # same csnaq_eval for every recorded run
python3 benchmarks/bench.py run --suite fit --label my-idea --extra-args newton   # fit suite: optimizer
python3 benchmarks/bench.py run --suite standard --cases n10_h1,n20_h2   # a subset of cases
```

## What is measured

The suites are defined in `suites.json`. In the `quick` and `standard` suites, every case uses about 1000
quartets. `propQuartets` and `propQuartetsFinal` are both set to `min(1, 1000 / binomial(ntaxa, 4))`
(the suite's `nquartets`), so the 10-taxon case uses all 210 of its quartets. A case can override
this with its own `propQuartets` and `propQuartetsFinal`.


| Suite | What it runs |
|---|---|
| `quick` | Three `snaq` cases, from 10 to 50 taxa, with the settings of `standard`. |
| `standard` | Eight `snaq` cases: 10 to 50 taxa, 1 to 3 reticulations, 100 runs and `Nfail` 40 each. Cases with 20 or more taxa use `--lazy`. |
| `fit` | Parameter optimization of the true network (`fitnumericalparameters!`, maxeval 1000), on 10 to 50 taxa. |

For every case, the history records:

| Metric | Meaning |
|---|---|
| `time_s` | C++: wall time of the whole `snaq` process. SNaQ.jl: time to build the `DataCF` plus the `snaq!` call, after a warm-up call, so Julia's start-up and compilation are excluded. `process_wall_s` is the whole process for both. |
| `peak_rss_mb` | Peak resident set size of the process. For SNaQ.jl this includes the Julia runtime; `base_rss_mb` is the size after loading the packages. |
| `best_score`, `run_scores`, `median_run_score` | As reported in the `.out` file. |
| `hwcd` | Hardwired-cluster distance between the best network and the true network (0 means the same topology). |
| `score_all_quartets` | SNaQ score of the best network, recomputed on the observed CFs of all quartets by `csnaq_eval`. This makes runs that used `propQuartets` < 1 comparable. |
| `rescored` | `hwcd`, `score_all_quartets` and `score_truth_all_quartets` recomputed by `bench.py rescore` from the saved `.out` file. The report uses these when present, so that all versions are evaluated by the same `csnaq_eval`. Run `rescore` after a change to `csnaq_eval` or to the expected CFs. |
| `score_refit_all_quartets` | SNaQ score of the best network's topology after fitting its parameters to all quartets. Compared with `score_truth_all_quartets`, it measures the topology alone. |
| `score_truth_all_quartets` | SNaQ score of the true network, after fitting its parameters to all quartets. |
| `fit_s`, `equations_s`, `score` | Fit suite: best of 3 repeats. |

Each record also holds:
- the label, implementation, git commit and whether the tree had uncommitted changes;
- the host (OS, CPU, cores), the thread count, the extra arguments and free-text `--notes`.

## Files

| Path | Content |
|---|---|
| `bench.py` | The framework: run suites, list runs, write reports. |
| `suites.json` | Suite and case definitions. Add cases here. |
| `datasets/<name>/` | `truth.net`, `genetrees.tre`, `start.tre` for each dataset. They are generated once by `julia/make_datasets.jl` and committed, so every run uses identical inputs. |
| `julia/` | The SNaQ.jl runners and the dataset generator. |
| `history/results.jsonl` | One JSON record per case per run. **Commit it:** it is the history. |
| `history/outputs/<run id>/` | The `.out` and `.networks` files of each case. |
| `reports/` | `REPORT.md` and SVG figures, written by `bench.py report`. Hovering over a bar or point in a browser shows its value. |

## Adding a case or a dataset

- **New case on an existing dataset:** add an object to a suite in `suites.json`. The
  available keys are `name`, `dataset`, `hmax`, `runs`, `Nfail`, `seed`, `propQuartets`,
  `propQuartetsFinal` and `lazy`. A suite-level `nquartets` sets both proportions for every case
  that does not set them (`maxeval` for the fit suite).
- **New dataset:** add it to `DATASETS` in `julia/make_datasets.jl` and run
  `julia --project=tools/julia benchmarks/julia/make_datasets.jl`.

Changing a case changes what its results mean. Give the changed case a new name, so the
history stays comparable.
