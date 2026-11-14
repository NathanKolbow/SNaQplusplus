# CSNaQ

C++ port of the `snaq!` function of [SNaQ.jl](https://github.com/JuliaPhylo/SNaQ.jl)
(branch `treechild-galled`, pull request #38): inference of phylogenetic networks
from gene trees by maximum pseudolikelihood using quartet concordance factors.
It also ports the parts of [PhyloNetworks.jl](https://github.com/JuliaPhylo/PhyloNetworks.jl)
that `snaq!` uses (extended Newick reading/writing and network manipulation).

The build produces one program, `snaq`.

## Requirements

- CMake ≥ 3.16
- A C++17 compiler with `<filesystem>` support: GCC ≥ 9, Clang ≥ 9, Apple Clang on macOS ≥ 10.15, or MSVC 2019+
- Git and an internet connection during the first configure step: CMake downloads and builds
  [NLopt](https://github.com/stevengj/nlopt) (v2.10.0) automatically. Nothing needs to be
  installed beforehand.

## Installation

Download the code, then build it with one command.

```sh
git clone https://github.com/NathanKolbow/CSNaQ.git
cd CSNaQ
make                      # macOS / Linux
```

On Windows, or on any system without `make`, run this instead from the `CSNaQ` folder:

```sh
cmake -P build.cmake
```

Both put the program in `bin/` (`bin/snaq`, or `bin/snaq.exe` on Windows).
Running `bin/snaq` with no arguments prints a short usage message.

Other `make` targets:

| Command | Effect |
|---|---|
| `make install` | install `snaq` to `~/.local/bin` (set `PREFIX=...` to change) |
| `make test` | build and run the unit tests |
| `make test-julia` | the unit tests, plus the tests that check `snaq`'s output files with SNaQ.jl |
| `make clean` | remove `build/` and `bin/` |

The same build with CMake directly:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
cmake --install build --prefix ~/.local    # optional
```

Optional CMake settings:

| Option | Default | Effect |
|---|---|---|
| `-DCSNAQ_NATIVE=ON` | `OFF` | Compile with `-march=native` |
| `-DCSNAQ_BUILD_TESTS=ON` | `OFF` | Build the unit tests (downloads GoogleTest) |
| `-DCSNAQ_JULIA_TESTS=ON` | `OFF` | With the unit tests, also run the tests that check `snaq`'s output files with SNaQ.jl (needs `julia` and the environment in `tools/julia`) |
| `-DCSNAQ_BUILD_TOOLS=ON` | `OFF` | Build the validation and benchmark tools (`csnaq_check`, `csnaq_classcheck`, `csnaq_bench_fit`, `csnaq_eval`) |

## Usage

```sh
snaq -i genetrees.tre -o results/net1 --start-tree start.tre --hmax 1 --runs 10 --Nfail 40
```

- `-i` / `--input`: gene trees, one Newick tree per line.
- `--start-tree`: starting topology. Either a Newick string (`"((a,b),c,(d,e));"`) or a file
  with one or more Newick lines. One topology is used by every run. Several topologies are
  used one per run, so give exactly one or `--runs` of them.
- `-o` / `--filename`: root name of the output files (the `filename` argument of `snaq!`).
  Default `snaq`. With `-o ""`, no files are written.

Every keyword argument of `snaq!` is available as a flag with the same name and default
value. The one exception is `ρ`, which is `--rho`. Run `snaq --help` for the full list, including:

| Flag | `snaq!` keyword | Default |
|---|---|---|
| `--hmax` | `hmax` | 1 |
| `--runs` | `runs` | 100 |
| `--Nfail` | `Nfail` | 50 |
| `--outgroup` | `outgroup` | none |
| `--seed` | `seed` | random |
| `--propQuartets`, `--propQuartetsFinal` | same | 1.0 |
| `--probST`, `--probQR` | same | 0.9, 0.0 |
| `--updateBL`, `--qinfTest` | same | true, false |
| `--ftolRel`, `--ftolAbs`, `--xtolRel`, `--xtolAbs` | same | 1e-8 |
| `--rho` | `ρ` | 0.0 |

- **Boolean flags:** given alone they mean `true`. They also accept a value
  (`--updateBL false` or `--updateBL=false`).
- **Restrictions:** these flags mirror `restrictionset`: `--max_level`, `--galled_tree`,
  `--galled_network`, `--rooted_tree_child`, `--weakly_tree_child` and `--strongly_tree_child`.
- **Advanced:** these are search keyword arguments that `snaq!` passes through:
  `--maxeval`, `--liktolAbs`, `--liktolRel`, `--opt_maxeval` and `--logprefix`.
- **Not in SNaQ.jl:** these flags control what `snaq` does differently from SNaQ.jl (see
  [Differences from SNaQ.jl](#differences-from-snaqjl)). The values in the last column
  reproduce SNaQ.jl's behavior.

  | Flag | Default | As SNaQ.jl |
  |---|---|---|
  | `--optimizer` | `newton` | `lbfgs` |
  | `--boundPrune` | true | false |
  | `--dedup` | true | false |
  | `--dedupCounts` | false | |
  | `--guideHybrids` | true | false |
  | `--minAddHybrid` | `Nfail` | 0 |
  | `--sweep` | false | false |
  | `--geneWeights` | false | false |
- **Other flags:**
  - `--lazy` computes observed CFs on demand, like `DataCF(trees; lazy=true)`. It requires
    `--propQuartets` < 1. It is on by default when `--propQuartets` and `--propQuartetsFinal`
    are both below 1, and off otherwise; `--lazy false` turns it off. The default stays off
    with `--qinfTest`, or when the starting topology does not have exactly the taxa of the
    gene trees, because the lazy `DataCF` supports neither. As in SNaQ.jl, the two differ for
    a quartet that a gene tree leaves unresolved: `--lazy` counts it for one of its tied
    topologies, and without `--lazy` it is skipped.
  - `--threads N` sets the number of threads; the default is all cores. Independent runs
    execute in parallel.

## Output

The files have the same names and layout as those written by `snaq!`:

| File | Content |
|---|---|
| `<filename>.out` | best network and its SNaQ score, then the networks of all runs |
| `<filename>.networks` | networks of all runs, sorted by SNaQ score |
| `<filename>.log` | settings and summary of the analysis |
| `<filename>_runs/run<j>.log` | log of run `j` |

The best network and its score are also printed to standard output.

When `--propQuartetsFinal` is between 0 and 1, or with `--lazy`, the scores in these files
are the scores of all networks on one sample of the quartets shared by all runs (a
proportion `--propQuartetsFinal` of them). Each run otherwise scores its network on its own
sample, and such scores are not comparable.

## Differences from SNaQ.jl

- **Random numbers.** The random number generator differs, so the same seed gives different
  runs in SNaQ.jl and in `snaq`.
- **Parallelism.** Independent runs execute in parallel threads.
- **Output files with `--lazy`.** With a lazy `DataCF`, `snaq!` returns the network that
  scores best on the quartets shared by all runs, but writes `.out` and `.networks` sorted by
  the score each network got in its own run, on that run's own sample of quartets. `snaq`
  writes both files sorted by the scores on the shared quartets, so their first network is
  the selected one.
- **Copied quirks.** Some behaviors of SNaQ.jl are copied on purpose. Each is marked in the
  source with a `// QUIRK:` comment.
- **Quartet classes.** Quartets whose four lineages enter the network through the same four
  cut edges have the same expected CFs. `snaq` builds one equation per such class and fits to
  the sum of the class's observed CFs. For each proposed network, the equations of the classes
  that the move leaves untouched are taken from the current network, for every kind of move
  (SNaQ.jl reuses equations after rNNI moves only).
- **Corrected expected CFs.** SNaQ.jl's equations give wrong expected CFs for some quartets of
  some networks: when a reticulation that cannot affect a quartet is removed and this leaves a
  root of degree 2, the root's two edges are fused and one of their two lengths is dropped from
  the equation. Class equations do not have this error. `csnaq_classcheck` compares both
  against a brute-force calculation: on `checks.net`, 88 of 133,552 tested quartets are wrong
  with SNaQ.jl's equations and none with class equations, and the scores from the two sets of
  equations differ for 64 of the 669 networks (on up to 3000 quartets per network).
- **Equations of reticulate classes.** SNaQ.jl builds the equation of a quartet that
  reticulations affect by copying the quartet's network once per way its lineages can pass
  the lowest hybrid node, editing each copy, and recursing. `snaq` leaves the network
  unchanged and moves the four lineages up it (`src/porteqn.hpp`). The equations have the same
  value for all parameter values, and are never larger: on `checks.net`,
  `csnaq_classcheck` finds no difference in score or gradient between the two builders
  (largest relative difference 1e-11), and 6% fewer equation nodes in total.
- **Evaluation of reticulate classes.** SNaQ.jl evaluates the equation of each quartet by
  recursion, with a list of the parameters seen on the way down from which it obtains the
  gradient. `snaq` stores every distinct sub-equation of a network's reticulate classes once
  (`src/eqndag.hpp`) and computes its expected CFs and their derivatives once per evaluation;
  classes with the same equation are evaluated together. The values are the same up to
  rounding: `csnaq_classcheck` compares the score, gradient and curvature with those of one
  recursion per class.
- **Networks are compared by the distance of their score from 0.** The search accepts a
  proposal, and `snaq` picks the best run, by how close the score is to 0 instead of by how
  high it is. The score is at most 0 up to rounding (expected CFs are floored at 1e-9), so
  both orders are the same unless a score is positive.
- **Parameter optimization.** By default, branch lengths and γ are optimized with a projected
  Newton method run to convergence: until a step close to a full Newton step is predicted to
  gain less than the tolerances. For networks, L-BFGS continues from its result in the rare
  fits where it stops for another reason (no step improves the score, or several steps in a
  row gain less than the tolerances).
  Before a fit stops, branches longer than 5 are checked separately: the derivative of the
  score with respect to such a length is close to 0 whatever a shorter branch would gain, so
  they are moved by a Newton step in exp(-length) when that step is predicted to raise the
  score, and the fit continues from there if it does.
  SNaQ.jl uses L-BFGS limited to `max(30, ntaxa)` evaluations; `--optimizer lbfgs` does the same.
- **Proposals rejected without a fit** (`--boundPrune`). An upper bound of the score of a
  proposed network is computed first: each class is given the parameter values that suit it
  best. If the bound cannot be accepted, the network is rejected without being fitted. It
  still counts as a failed proposal, so this does not change the result of a search.
- **Topologies are not fitted twice** (`--dedup`). A run recognizes the topologies it has
  already fitted. Such a proposal is not fitted again, and by default not counted as a failed
  proposal either: another proposal is drawn, so a run stops after `Nfail` distinct
  topologies failed in a row (`--dedupCounts` counts them instead). When `--propQuartets` is
  1, a topology fitted by one run is not fitted again by the others.
- **Hybrid additions.** SNaQ.jl proposes a hybrid addition 5% of the time, between two edges
  drawn uniformly. In `snaq`:
  - half of them are drawn between edges where the observed CFs exceed the current network's
    expected CFs (`--guideHybrids`);
  - a run whose network has fewer than `hmax` reticulations does not stop before `Nfail`
    hybrid additions were evaluated since its last accepted move (`--minAddHybrid`).
- **Restrictions are checked on the network as it is written.** The tree-child conditions
  read the direction of the edges and whether the root can be placed on them. SNaQ.jl checks
  a proposed network as the move left it, where these are not up to date for every edge;
  ported as is, this accepted networks that do not meet `--weakly_tree_child` once they are
  read from the output files. `snaq` computes them again before the check.
- **Sweep of the local moves** (`--sweep`, off by default). A run of SNaQ.jl stops after
  `Nfail` random proposals failed in a row, whether or not a move that it did not draw would
  improve its network. With `--sweep`, a run that would stop first proposes every local move
  of its network: the NNIs and, for each reticulation, the local moves of its origin and of
  its target and the change of its direction. If one is accepted, the run goes on with random
  proposals; it stops when none is.
- **Gene weights** (`--geneWeights`, off by default). Each quartet is weighted by the
  proportion of gene trees that its observed CFs come from. SNaQ.jl weighs all quartets
  equally. The two are the same when every gene tree has every taxon.
- **Input.** Gene trees must be in Newick format. CF tables, NEXUS files and multiple
  alleles are not supported.
- **Tests.** The tests that were ported, and those that were not, are listed in
  `TESTS_SKIPPED.md`.

## Validation tools

| Path | Purpose |
|---|---|
| `tools/julia/setup.jl` | sets up a Julia environment with SNaQ.jl, expected at `../SNaQ.jl` relative to the repository root |
| `tools/csnaq_check.cpp`, `tools/julia/check_networks.jl` | compare expected CFs, SNaQ scores and gradients of the per-quartet equations with SNaQ.jl on the networks of a Newick file (e.g. `checks.net`) |
| `tools/classcheck.cpp`, `tools/bruteforce_cf.hpp` | compare the quartet-class likelihood with the per-quartet one, both with brute-force expected CFs, the two builders of reticulate class equations, and the evaluation on shared sub-equations with one recursion per class: `build/csnaq_classcheck checks.net 3000` |
| `validation/checks_net_comparison.csv` | result of that comparison on `checks.net` |
| `tools/julia/check_restrictions.jl` | runs `snaq` with each restriction flag and checks every network it wrote with SNaQ.jl's `restrictionset` (on a copy, after `semidirectnetwork!`); run by `make test-julia` |
| `tools/bench_fit.cpp`, `tools/julia/bench_fit.jl` | time numerical parameter optimization in both implementations, with peak memory |
| `tools/julia/make_test_fixtures.jl` | regenerates the test fixtures in `tests/data/` |

## License

MIT; see [LICENSE](LICENSE). This includes the license notices of SNaQ.jl and PhyloNetworks.jl.
NLopt is downloaded at build time and linked statically. Its L-BFGS code is licensed under
the GNU LGPL v2.1 or later.
