# Ported and skipped tests

This file lists every test file of SNaQ.jl (branch `treechild-galled`) and how it was ported
to the C++ test suite in `tests/` (GoogleTest; run with `make test`). Only functionality
needed by `snaq!` was ported (see the README), so tests of other SNaQ.jl functions are skipped.

## General notes

- **Random inputs.** SNaQ.jl's tests draw them with Julia's RNG and PhyloCoalSimulations.jl
  (simulated gene trees, `generate_net`, `generate_tree`). Those inputs are generated once by
  `tools/julia/make_test_fixtures.jl` and stored in `tests/data/`.
- **Moves and optimizer perturbations.** Random choices made during moves or perturbations
  use CSNaQ's own RNG, so sequences of random moves differ from Julia's.
- **Warnings.** Assertions on Julia warnings (`@test_logs`, `@test_nowarn`) are checked as
  "does not throw".

## Test files

| SNaQ.jl test file | C++ test file | Status |
|---|---|---|
| `test_gradient_opt/test_CF_recursive_blocks.jl` | `test_cf_recursive_blocks.cpp` | Ported. Not ported: `coalmask` (an internal helper CSNaQ does not use) and direct calls to the internal recursion `computeexpectedCFandgradientrecur!`. Its values and derivatives are checked through `computeexpectedCF` and finite differences instead. |
| `test_gradient_opt/test_opt_API.jl` | `test_opt_api.cpp` | Ported. |
| `test_gradient_opt/test_search_API.jl` | `test_search_api.cpp` | Ported. |
| `test_inplace_updates/test_rNNI_inplace.jl` | `test_inplace_and_data.cpp` | Ported with changes; see below. |
| `test_perfectData.jl` | `test_inplace_and_data.cpp` | Ported. |
| `test_negativeEdgeLengths.jl` | `test_inplace_and_data.cpp` | Ported. |
| `test_propQuartets.jl` | `test_propquartets_lazy_bugs.cpp` | Ported. Not ported: assertions on the lazy `DataCF` object itself (`computeSNaQscore!(truetre, ldcf)` erroring when no CFs are computed, and `length(ldcf.quartet)`), because the C++ port has no `DataCF` object. |
| `test_lazy_datacf.jl` | `test_propquartets_lazy_bugs.cpp` | Partly ported: lazy and non-lazy observed CFs give the same SNaQ score, the taxa and number of quartets, the quartet of rank 1, and an error past the last quartet. Not ported: writing and reading lazy CFs to and from CSV (`write`, `LazyDataCF(file)`), `computeSNaQscore!` on a lazy `DataCF` with `numQuartets`, display methods, and the safeguards of functions that do not exist in the port. |
| `test_bugs.jl` | `test_propquartets_lazy_bugs.cpp` | Partly ported: "log file describes the restrictions used" (the BUCKy CF table is converted to a fixture by Julia) and "hmax enforced on starting topologies". Not ported: `loglik`/`loglik!` (SNaQ v1 scale), expected CFs stored in a `DataCF`, `fittedquartetCF`, `compositedeviance`/`compositeloglik`, and reading v1/v2 output files (`readsnaqnetwork`, `readallsnaqnetworks`). |
| `test_network_moves/test_add_remove_retic.jl` | `test_network_moves.cpp` | Ported. The move SNaQ.jl samples with `Random.seed!(0)` is hard-coded and gives the same network as in Julia. |
| `test_network_moves/test_move_target_origin.jl` | `test_network_moves.cpp` | Ported. |
| `test_network_moves/test_rSPR_moves.jl` | `test_network_moves.cpp` | Ported. |
| `test_network_moves/test_rNNI_moves.jl` | `test_network_moves.cpp` | Ported. `performrandomrNNI!` (a test-only helper) is reimplemented in the test. |
| `test_network_moves/test_all_moves.jl` | `test_network_moves.cpp` | Ported. `generate_tree` is replaced by a random tree generator in `tests/test_support.cpp`. |
| `test_confint_qcf_genetrees.jl` | none | Skipped: `confintqCF_genetrees` is not part of the port. |
| `test_bootstrap.jl` | none | Skipped: `bootsnaq` is not part of the port. |
| `test_multipleAlleles.jl` | none | Skipped: multiple alleles and CF tables are not part of the port. |
| `test_readInputData.jl` | none | Skipped: `readtrees2CF`, `readtableCF`, `tablequartetCF` and quartet files are not part of the port. Observed CFs from gene trees are tested by `LazyDataCF.MatchesNonLazy` and `PropQuartets.LazyDataCF`. |
| Aqua.jl code-quality checks (`runtests.jl`) | none | Not applicable to C++. |

## `test_rNNI_inplace.jl`

This file is not part of SNaQ.jl's `runtests.jl`, and it does not run in SNaQ.jl:

- **API error.** It calls `computeSNaQscore!` with an outdated argument type.
- **KeyError when fixed.** With that type fixed, it raises a `KeyError`, because it passes
  the parameter map of the network *before* the move.

Ported with three changes:

1. **Parameter map.** The test passes the map of the network *after* the move, as SNaQ.jl's
   `search` does.
2. **Changed parameter layout.** Moves that change which edges carry parameters are not
   compared. An example is a pendant edge becoming internal while the parameter count stays
   the same. CSNaQ's search rebuilds the equations for such moves.
3. **rNNI(1) with `v` the parent of `u`.** These moves are not compared. For them, SNaQ.jl's
   in-place rNNI(1) update gives equations that differ from a full rebuild, because
   quartets with taxa below `t` are not re-derived. CSNaQ's search rebuilds them.

Changes 2 and 3 concern SNaQ.jl's in-place updates, which are ported and tested
(`applyrNNI1update`, `applyrNNI2update`). CSNaQ's search no longer uses them: it works on
quartet classes, and reuses the equations of the classes a move leaves untouched (see the
README).

## Tests that are not ported from SNaQ.jl

| File | What it tests |
|---|---|
| `tests/test_classes.cpp` | Quartet classes: expected CFs against a brute-force calculation, the score against the per-quartet sum, gradients against finite differences, and equations reused after random moves against equations built from scratch. |
| `tests/test_search_shortcuts.cpp` | The score upper bound (never below a fitted score; does not change a search), the topology signature, the Newton optimizer against L-BFGS run to convergence, the candidate pairs for hybrid additions, and gene weights. |

The tests ported from SNaQ.jl compute expected CFs and scores of networks with
`computeexpectedCFmatrix` and `computeSNaQscore_net`, which now use quartet classes. For
networks where SNaQ.jl's per-quartet equations are wrong (see "Corrected expected CFs" in the
README), these values differ from SNaQ.jl's.
