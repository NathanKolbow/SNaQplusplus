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
| `-DCSNAQ_BUILD_TOOLS=ON` | `OFF` | Build the validation and benchmark tools (`csnaq_check`, `csnaq_bench_fit`) |

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
- **Other flags:**
  - `--lazy` computes observed CFs on demand, like `DataCF(trees; lazy=true)`. It requires
    `--propQuartets` < 1.
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

## Differences from SNaQ.jl

- **Random numbers.** The random number generator differs, so the same seed gives different
  runs in SNaQ.jl and in `snaq`.
- **Parallelism.** Independent runs execute in parallel threads.
- **Copied quirks.** Some behaviors of SNaQ.jl are copied on purpose. Each is marked in the
  source with a `// QUIRK:` comment.
- **Input.** Gene trees must be in Newick format. CF tables, NEXUS files and multiple
  alleles are not supported.

## Validation tools

| Path | Purpose |
|---|---|
| `tools/julia/setup.jl` | sets up a Julia environment with SNaQ.jl, expected at `../SNaQ.jl` relative to the repository root |
| `tools/csnaq_check.cpp`, `tools/julia/check_networks.jl` | compare expected CFs, SNaQ scores and gradients with SNaQ.jl on the networks of a Newick file (e.g. `checks.net`) |
| `validation/checks_net_comparison.csv` | result of that comparison on `checks.net` |
| `tools/bench_fit.cpp`, `tools/julia/bench_fit.jl` | time numerical parameter optimization in both implementations |

## License

MIT; see [LICENSE](LICENSE). This includes the license notices of SNaQ.jl and PhyloNetworks.jl.
NLopt is downloaded at build time and linked statically. Its L-BFGS code is licensed under
the GNU LGPL v2.1 or later.
