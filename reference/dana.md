# Differential analysis (dana)

Feature-wise [`stats::lm()`](https://rdrr.io/r/stats/lm.html) or
[`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html) models of an
omics data matrix. Supports likelihood ratio tests (LRT) and parallel
computation.

## Usage

``` r
dana(
  X,
  sample_data,
  formula_rhs,
  term_LRT = NULL,
  model_control = list(),
  platform = c("ms", "nmr", "ngs"),
  assay = NULL,
  verbose = TRUE
)
```

## Arguments

- X:

  A numeric matrix with samples in rows and features in columns. Sample
  IDs in row names must match the format from `sample_id` column in
  `sample_data`.

- sample_data:

  A data frame containing sample-level data. Must have a `sample_id`
  column matching row names in `X` and `sample_data`.

- formula_rhs:

  A one-sided formula (e.g., `~ group + (1|subject)`). Must not contain
  a response variable.

- term_LRT:

  Optional. Character vector of formula terms to test via LRT. Random
  effects must be written without parentheses (e.g., `"1 | group"`).

- model_control:

  Optional. List of control arguments passed to the model.

- platform:

  Character string indicating the omics platform (e.g., `"ms"`, `"nmr"`,
  `"ngs"`).

- assay:

  Optional. Character string indicating the name of the platform assay
  (e.g., `"lipidomics"`).

- verbose:

  Logical. If TRUE, prints progress messages.

## Value

An object of class `"dana"`:

- X:

  Matched data matrix.

- sdata:

  Matched sample data.

- fit:

  Data frame of model coefficients and confidence intervals per feature.

- lrt:

  Likelihood ratio test results (if `term_LRT` is specified).

- ranef:

  Random effects variance components (if using mixed models).

- errors:

  A data frame logging any model fitting errors per feature.

## Details

Models are fit independently for each feature using
[`stats::lm()`](https://rdrr.io/r/stats/lm.html) or
[`lmerTest::lmer()`](https://rdrr.io/pkg/lmerTest/man/lmer.html),
depending on whether `dana()` detects random effects in `formula_rhs`.
Feature-wise models can be evaluated in parallel using
[`future::plan()`](https://future.futureverse.org/reference/plan.html),
with optional progress updates via
[`progressr::with_progress()`](https://progressr.futureverse.org/reference/with_progress.html).

## See also

[`stats::lm()`](https://rdrr.io/r/stats/lm.html),
[`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html),
[`lmerTest::lmer()`](https://rdrr.io/pkg/lmerTest/man/lmer.html)
parameters.

## Examples

``` r
mock_X <- matrix(
  rnorm(50 * 10) +
    rep(c(rep(0, 25), rep(2, 25)), each = 10) * rep(1:10 %in% 1:3, each = 50),
  nrow = 50
)

rownames(mock_X) <- paste0("sample", 1:50)
colnames(mock_X) <- paste0("feat", 1:10)

sample_data <- data.frame(
  sample_id = rownames(mock_X),
  group = factor(rep(c("A", "B"), each = 25)),
  subject = factor(rep(1:25, each = 2)),
  row.names = rownames(mock_X)
)

# Example with parallel computation setup (not run)
# future::plan(multisession)
# progressr::handlers(global = TRUE)
# progressr::with_progress({
  result <- dana(X = mock_X,
                 sample_data = sample_data,
                 formula_rhs = ~ group + (1 | subject),
                 term_LRT = c("group", "1 | subject"), # Multiple terms allowed
                 platform = "ms",
                 assay = "lipidomics",
                 verbose = FALSE
                 )
#> Computing profile confidence intervals ...
#> refitting model(s) with ML (instead of REML)
#> refitting model(s) with ML (instead of REML)
#> boundary (singular) fit: see help('isSingular')
#> Computing profile confidence intervals ...
#> boundary (singular) fit: see help('isSingular')
#> refitting model(s) with ML (instead of REML)
#> Warning: convergence code 3 from bobyqa: bobyqa -- a trust region step failed to reduce q
#> refitting model(s) with ML (instead of REML)
#> Computing profile confidence intervals ...
#> refitting model(s) with ML (instead of REML)
#> refitting model(s) with ML (instead of REML)
#> Computing profile confidence intervals ...
#> refitting model(s) with ML (instead of REML)
#> refitting model(s) with ML (instead of REML)
#> Computing profile confidence intervals ...
#> refitting model(s) with ML (instead of REML)
#> refitting model(s) with ML (instead of REML)
#> boundary (singular) fit: see help('isSingular')
#> Computing profile confidence intervals ...
#> boundary (singular) fit: see help('isSingular')
#> refitting model(s) with ML (instead of REML)
#> refitting model(s) with ML (instead of REML)
#> boundary (singular) fit: see help('isSingular')
#> Computing profile confidence intervals ...
#> boundary (singular) fit: see help('isSingular')
#> refitting model(s) with ML (instead of REML)
#> refitting model(s) with ML (instead of REML)
#> Computing profile confidence intervals ...
#> refitting model(s) with ML (instead of REML)
#> refitting model(s) with ML (instead of REML)
#> boundary (singular) fit: see help('isSingular')
#> Computing profile confidence intervals ...
#> Warning: convergence code 3 from bobyqa: bobyqa -- a trust region step failed to reduce q
#> boundary (singular) fit: see help('isSingular')
#> refitting model(s) with ML (instead of REML)
#> Warning: convergence code 3 from bobyqa: bobyqa -- a trust region step failed to reduce q
#> refitting model(s) with ML (instead of REML)
#> Warning: convergence code 3 from bobyqa: bobyqa -- a trust region step failed to reduce q
#> Computing profile confidence intervals ...
#> refitting model(s) with ML (instead of REML)
#> refitting model(s) with ML (instead of REML)
# })

# Modify `dana` object at once with pipes (not run)
# dana_obj <- dana_obj |> adjust_pval() |> add_feat_name() |> ready_plots()
```
