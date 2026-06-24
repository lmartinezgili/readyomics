# PERMANOVA with flexible permutation control

Performs PERMANOVA (Permutational Multivariate Analysis of Variance).
Supports both joint-term (default
[`vegan::adonis2()`](https://vegandevs.github.io/vegan/reference/adonis.html))
and single-term testing when `independent = TRUE`. Several distance
methods, and fine-grained permutation control.

## Usage

``` r
permanova(
  X,
  sample_data,
  formula_rhs,
  dist_control = list(method = "euclidean", diag = FALSE, upper = FALSE),
  perm_control = list(joint_terms = list(control = permute::how(blocks = NULL, nperm =
    999))),
  independent = TRUE,
  platform = c("ms", "nmr", "ngs"),
  assay = NULL,
  seed = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- X:

  A processed matrix or data frame of features (samples in rows,
  features in columns).

- sample_data:

  A `data.frame` containing sample-level data. Row names must match
  those in `X`.

- formula_rhs:

  A one-sided formula (e.g., `~ group + age`).

- dist_control:

  A named list of arguments to control distance calculation. Must
  contain at least `method`. Defaults to `"Euclidean"` via
  [`stats::dist()`](https://rdrr.io/r/stats/dist.html).

- perm_control:

  A named list specifying
  [`permute::shuffleSet()`](https://rdrr.io/pkg/permute/man/shuffleSet.html)
  parameters. By default, `joint_terms` parameters will be used, with
  same
  [`vegan::adonis2()`](https://vegandevs.github.io/vegan/reference/adonis.html)
  defaults, unless variable-specific permutation settings are added as
  named list elements (e.g. perm_control = list(joint_terms = , age = ,
  sex = )).

- independent:

  Logical. If `TRUE`, a PERMANOVA test for each variable in
  `formula_rhs` is performed.

- platform:

  A string specifying the omics platform (`"ms"`, `"nmr"`, `"ngs"`).
  Used for annotation.

- assay:

  Optional. Character string giving the assay name for annotation (e.g.,
  `"lipidomics"`).

- seed:

  Optional integer. If provided, sets the random seed for reproducible
  permutation results.

- verbose:

  Logical. If `TRUE`, prints diagnostic messages.

- ...:

  Additional arguments passed to
  [`vegan::adonis2()`](https://vegandevs.github.io/vegan/reference/adonis.html).

## Value

A named `list` with three elements:

- X_dist:

  A `dist` object.

- perm_matrix_joint:

  A `matrix` from
  [`permute::shuffleSet()`](https://rdrr.io/pkg/permute/man/shuffleSet.html)
  joint_terms control.

- permanova_joint:

  A `data.frame` of PERMANOVA results using the full model.

- permanova_indep:

  A `data.frame` of a PERMANOVA results for each predictor, or `NULL` if
  `independent = FALSE`.

## Details

- Supports both [`stats::dist()`](https://rdrr.io/r/stats/dist.html) and
  [`vegan::vegdist()`](https://vegandevs.github.io/vegan/reference/vegdist.html)
  for distance matrix computation.

- Distance method must be specified in `dist_control$method`.

- Permutation design is controlled via the permute package using
  [`permute::shuffleSet()`](https://rdrr.io/pkg/permute/man/shuffleSet.html).

- If `seed` is supplied, the same permutations will be used across runs
  for reproducibility.

## See also

- [`stats::dist()`](https://rdrr.io/r/stats/dist.html) and
  [`vegan::vegdist()`](https://vegandevs.github.io/vegan/reference/vegdist.html)
  for information on available distances.

- [`vegan::adonis2()`](https://vegandevs.github.io/vegan/reference/adonis.html)
  and
  [`permute::shuffleSet()`](https://rdrr.io/pkg/permute/man/shuffleSet.html)
  for control options and details.

- [`process_ngs()`](https://lmartinezgili.github.io/readyomics/reference/process_ngs.md)
  to pre-process and normalize an `X` NGS dataset.

- [`process_ms()`](https://lmartinezgili.github.io/readyomics/reference/process_ms.md)
  to pre-process and normalize an `X` MS dataset.

## Examples

``` r
# Mock data
X <- matrix(rnorm(40), nrow = 10,
            dimnames = list(paste0("sample", 1:10),
                            paste0("feat", 1:4)))
sample_data <- data.frame(
  sample_id = rownames(X),
  group = factor(rep(c("A", "B"), each = 5)),
  age = rep(20:29, length.out = 10),
  row.names = rownames(X),
  stringsAsFactors = FALSE
)

# Simple control structures
dist_control <- list(method = "euclidean")
perm_control <- list(
  joint_terms = list(control = permute::how(blocks = NULL, nperm = 9)),
  group = list(control = permute::how(blocks = NULL, nperm = 9)),
  age = list(control = permute::how(blocks = NULL, nperm = 9))
)

result <- permanova(
  X = X,
  sample_data = sample_data,
  formula_rhs = ~ group + age,
  dist_control = dist_control,
  perm_control = perm_control,
  independent = TRUE,
  platform = "ms",
  assay = "lipidomics",
  seed = 42,
  verbose = FALSE
)
```
