# treebalance 1.2.2

* **Changed results and return types:** The Wedderburn-Etherington numbers in `wedEth` are now
  stored exactly as big integers (`bigz`, package `gmp`). Previously they were stored as double,
  which made them inexact from n = 50 and infinite from n = 793 on. As a consequence,
  `furnasI()`, `furnasI_inv()`, `treenumber()`, `treenumber_inv()` and `getfurranks()` now
  compute exactly for all trees with up to 2545 leaves; previously they silently returned wrong
  results for n >= 51. `we_eth()` now returns `bigz` by default.
* **Changed results and return types:** `colPlaLab()` now computes the Colijn-Plazzotta rank
  exactly and returns it as `bigz` by default; previously the ranks were wrong from n = 10 and
  infinite from n = 14 on. `colPlaLab_inv()` and `cPL_inv()` now also compute exactly.
* New argument `type = "double"` in `we_eth()`, `furnasI()`, `treenumber()` and `colPlaLab()`
  to obtain the result as double. It is only available as long as all values can be
  represented exactly (n <= 48, n <= 47 and n <= 9, respectively); otherwise the functions stop.
* `furnasI_inv()`, `treenumber_inv()`, `colPlaLab_inv()` and `cPL_inv()` accept the rank as
  `bigz`, character or double (the latter up to 2^53).
* **Changed results:** `collesslikeI(dissim = "mdm")` now computes the mean deviation from the
  median as documented. Previously the sum instead of the mean was returned, i.e. the balance
  value of every vertex was multiplied by its number of children.
* `rQuartetI()` now takes the shape value `q0` into account (previously it was ignored, so results
  were only correct for `q0 = 0`) and checks that `shapeVal` has exactly 5 entries.
* `B1I()`, `totCophI()` and `areaPerPairI()` now return 0 (resp. the correct value) for star trees
  instead of `NA`.
* `collessI(method = "corrected")` now returns 0 for `n = 2` as documented instead of `NaN`.
* `is_binary()` now returns `FALSE` instead of `NULL` for non-binary trees, so the binary-only
  functions show their intended error message.
* `collesslikeI()` now also accepts functions for `f.size` and `dissim` directly, not only their
  names as strings (needed e.g. for functions defined within another function or package).
* `B2I()` and `sShapeI()` now stop for the invalid logarithm base 1 (and `sShapeI()` also for
  non-positive bases) instead of returning `Inf`/`NaN`.
* `tree_merge()` now also works for trees with edge lengths (previously it failed with
  "tree 'x' has no root edge") and keeps tip labels and edge lengths when two single leaves
  are merged.
* `furnasI_inv()` and `treenumber_inv()` now also work when `treebalance` is not attached,
  e.g. when used from within another package.
* Updated the citation to the published book "Tree Balance Indices -- A Comprehensive Survey"
  (Springer, 2023).

# treebalance 1.1.0

* Added a `NEWS.md` file to track changes to the package.
* Added arXiv identifier of "Tree balance indices: a comprehensive survey" to DESCRIPTION.
* Corrected a mistake in the function furnasI_inv().

# treebalance 1.2.0
* Added modified maximum difference in width.
* Changed FurnasI such that it can deal with larger numbers.
