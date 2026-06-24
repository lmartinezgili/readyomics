# Adjust P-values in a `dana` object

Applies multiple testing correction to P-values from differential
analysis results returned by the
[`dana()`](https://lmartinezgili.github.io/readyomics/reference/dana.md)
function. Supports multiple adjustment methods and both coefficient and
likelihood ratio test (LRT) P-values.

## Usage

``` r
adjust_pval(
  dana_obj,
  padj_by = c("all", "terms"),
  padj_method = NULL,
  padj_method_LRT = NULL,
  ihw_covar = NULL,
  ihw_covar_id = NULL,
  ihw_args = list(),
  storey_args = list(),
  verbose = TRUE
)
```

## Arguments

- dana_obj:

  A `dana` class object returned by the
  [`dana()`](https://lmartinezgili.github.io/readyomics/reference/dana.md)
  function.

- padj_by:

  Character string. Whether P-value adjustment should be done globally
  across all coefficients (`"all"`) or separately for each coefficient
  term (`"terms"`).

- padj_method:

  Character vector of one or more methods for adjusting P-values from
  coefficient tests. Defaults to `"BH"`.

- padj_method_LRT:

  Character vector of one or more methods for adjusting P-values from
  LRT tests. Defaults to `"BH"`. P-values from LRT tests will always be
  adjusted independently for each LRT term.

- ihw_covar:

  Data frame containing covariable(s) used for
  [`IHW::ihw()`](https://rdrr.io/pkg/IHW/man/ihw.default.html). Must
  contain a `"feat_id"` column, matching `dana` object `"feat_id"`
  labels.

- ihw_covar_id:

  Character string. Column name in `ihw_covar` used as `covariates`
  argument in
  [`IHW::ihw()`](https://rdrr.io/pkg/IHW/man/ihw.default.html).

- ihw_args:

  Named list. Additional arguments passed to
  [`IHW::ihw()`](https://rdrr.io/pkg/IHW/man/ihw.default.html). Do not
  provide `covariates` argument, as it will be obtained from
  `ihw_covar`.

- storey_args:

  Named list. Additional arguments passed to
  [`qvalue::qvalue()`](https://rdrr.io/pkg/qvalue/man/qvalue.html).

- verbose:

  Logical. Whether to print informative messages. Defaults to `TRUE`.

## Value

A modified `dana` object with new columns in the `$fit` and `$lrt` data
frames for each adjusted P-value method applied (e.g. `padj_BH`,
`padj_storey_group`).

## Details

Available adjustment methods include: `"BH"`, `"bonferroni"`, `"BY"`,
`"fdr"`, `"hochberg"`, `"holm"`, `"hommel"`, `"IHW"`, and `"storey"`.

## See also

- [`dana()`](https://lmartinezgili.github.io/readyomics/reference/dana.md)
  for fitting differential analysis models on omics datasets.

- [`IHW::ihw()`](https://rdrr.io/pkg/IHW/man/ihw.default.html) for
  inverted hypothesis weighting method details.

- [`qvalue::qvalue()`](https://rdrr.io/pkg/qvalue/man/qvalue.html) for
  Storey's qvalue method details.

## Examples

``` r
set.seed(123)
mock_X <- matrix(rnorm(20 * 5), nrow = 20)
colnames(mock_X) <- paste0("feat_", seq_len(5))
rownames(mock_X) <- paste0("sample_", seq_len(20))

sample_data <- data.frame(
  sample_id = rownames(mock_X),
  group = factor(rep(c("A", "B"), each = 10)),
  time = factor(rep(c("T1", "T2"), times = 10)),
  subject_id = factor(rep(seq_len(10), each = 2)),
  stringsAsFactors = FALSE
)
rownames(sample_data) <- sample_data$sample_id

fit_df <- data.frame(
  feat_id = rep(colnames(mock_X), each = 2),
  Coefficient = rep(c("(Intercept)", "groupB"), 5),
  Estimate = rnorm(10),
  `Pr(>|t|)` = runif(10),
  stringsAsFactors = FALSE
)

# Mock dana object
dana_obj <- list(
  X = mock_X,
  sdata = sample_data,
  formula_rhs = ~ group,
  fit = fit_df,
  lrt = data.frame(),
  ranef = data.frame()
)
class(dana_obj) <- "dana"

# Add adjusted P-values
dana_obj <- dana_obj |>
  adjust_pval(padj_method = c("BH", "bonferroni"),
              padj_method_LRT = NULL,
              padj_by = "terms",
              verbose = FALSE)
```
