# Process next generation sequencing data

This function performs quality control, filtering, normalization, and
transformation of sequencing data raw counts. It can also build
`phyloseq` objects for downstream ecological analyses, and optionally
returns intermediate processing steps.

## Usage

``` r
process_ngs(
  X,
  sample_data,
  taxa_table = NULL,
  phylo_tree = NULL,
  remove_ids = NULL,
  min_reads = 500,
  min_prev = 0.1,
  normalise = c("load", "TSS", "none"),
  load_colname = NULL,
  min_load = 10000,
  transform = c("clr", "log", "none"),
  impute_control = list(method = "GBM", output = "p-counts", z.delete = FALSE, z.warning
    = 1, suppress.print = TRUE),
  raw_phyloseq = TRUE,
  eco_phyloseq = TRUE,
  return_all = FALSE,
  verbose = TRUE
)
```

## Arguments

- X:

  A numeric matrix or data frame of raw counts with samples as rows and
  features (e.g., taxa) as columns. Row names must be sample IDs.

- sample_data:

  A data frame containing sample-level data. Must include a column named
  `sample_id` with matching row names with `X`.

- taxa_table:

  Optional. Taxonomy annotation table to build `phyloseq` objects. Row
  names must match column names of `X`.

- phylo_tree:

  Optional. Phylogenetic tree to add to `phyloseq` objects.

- remove_ids:

  A regex or character vector to filter rows in `X`. Set to `NULL` to
  skip.

- min_reads:

  Numeric. Minimum number of total reads required per sample. Default is
  500.

- min_prev:

  Numeric between 0 and 1. Minimum feature prevalence threshold. Default
  is 0.1 (i.e., feature must be present in \>= 10 % of samples).

- normalise:

  Normalization method. One of `"load"` (microbial load data), `"TSS"`
  (total sum scaling), or `"none"`.

- load_colname:

  Column name in `sample_data` containing microbial load values.
  Required if `normalise = "load"`.

- min_load:

  Numeric. Default is 1e4. Warns if any microbial load value \<
  min_load.

- transform:

  Transformation method. One of `"clr"` (centered log-ratio with zero
  imputation), `"log"` (pseudo-log using
  [`log1p()`](https://rdrr.io/r/base/Log.html)), or `"none"`. Note: When
  using `"clr"`, zero values are imputed using
  [`zCompositions::cmultRepl()`](https://rdrr.io/pkg/zCompositions/man/cmultRepl.html).

- impute_control:

  A named list of arguments to be passed to
  [`zCompositions::cmultRepl()`](https://rdrr.io/pkg/zCompositions/man/cmultRepl.html).

- raw_phyloseq:

  Logical. If `TRUE`, constructs a `phyloseq` object with the table of
  raw counts (filtered failed runs if needed). Default is `TRUE`.

- eco_phyloseq:

  Logical. If `TRUE`, constructs a `phyloseq` object with the ecosystem
  abundances (i.e. after `normalise = "load"`). Default is `TRUE`.

- return_all:

  Logical. If `TRUE`, additional intermediate data matrices
  (`X_matched`, `X_norm`, `X_prev`) are included in the output. Default
  is `FALSE`.

- verbose:

  Logical. If `TRUE`, prints progress messages during execution. Default
  is `TRUE`.

## Value

A named list containing:

- `X_processed`:

  Matrix of processed feature counts after filtering, normalization, and
  transformation.

- `sdata_final`:

  Matched and filtered `sample_data` corresponding to retained samples.

- `phyloseq_raw`:

  `phyloseq` object created from raw filtered data. `NULL` if
  `raw_phyloseq = FALSE`.

- `phyloseq_eco`:

  `phyloseq` object from ecosystem abundance data. `NULL` if
  `eco_phyloseq = FALSE` or `normalise != "load"`.

- `X_matched`:

  (Optional) Matched and filtered count matrix, pre-normalization.
  Returned only if `return_all = TRUE`.

- `X_norm`:

  (Optional) Normalized count matrix. Returned only if
  `return_all = TRUE`.

- `X_prev`:

  (Optional) Prevalence-filtered matrix, pre-transformation. Returned
  only if `return_all = TRUE`.

## Details

- Zeros are imputed with
  [`zCompositions::cmultRepl()`](https://rdrr.io/pkg/zCompositions/man/cmultRepl.html)
  before CLR transformation.

- QC or other samples are removed if `remove_ids` is specified.

- Sample IDs in `X` and `sample_data` row names are matched and aligned.

- Can generate both a `phyloseq_raw` phyloseq object containing raw
  counts and a `phyloseq_eco` object with ecosystem counts, if a
  `load_colname` column from `sample_data` is provided to normalize the
  counts by microbial load (recommended best practice).

## References

\#' McMurdie, P. J., & Holmes, S. (2013). phyloseq: An R package for
reproducible interactive analysis and graphics of microbiome census
data. *PLoS ONE*, 8(4), e61217.
[doi:10.1371/journal.pone.0061217](https://doi.org/10.1371/journal.pone.0061217)

Martín-Fernández, J. A., Hron, K., Templ, M., Filzmoser, P., &
Palarea-Albaladejo, J. (2015). Bayesian-multiplicative treatment of
count zeros in compositional data sets. *Statistical Modelling*, 15(2),
134–158.
[doi:10.1177/1471082X14535524](https://doi.org/10.1177/1471082X14535524)

Palarea-Albaladejo, J., & Martín-Fernández, J. A. (2015).
zCompositions—R package for multivariate imputation of left-censored
data under a compositional approach. *Chemometrics and Intelligent
Laboratory Systems*, 143, 85–96.
[doi:10.1016/j.chemolab.2015.02.019](https://doi.org/10.1016/j.chemolab.2015.02.019)

Gloor, G. B., Macklaim, J. M., Pawlowsky-Glahn, V., & Egozcue, J. J.
(2017). Microbiome datasets are compositional: And this is not optional.
*Frontiers in Microbiology*, 8, 2224.
[doi:10.3389/fmicb.2017.02224](https://doi.org/10.3389/fmicb.2017.02224)

Vandeputte, D., Kathagen, G., D’hoe, K., Vieira-Silva, S.,
Valles-Colomer, M., Sabino, J., Wang, J., Tito, R. Y., De Commer, L.,
Darzi, Y., Vermeire, S., Falony, G., & Raes, J. (2017). Quantitative
microbiome profiling links gut community variation to microbial load.
*Nature*, 551(7681), 507–511.
[doi:10.1038/nature24460](https://doi.org/10.1038/nature24460)

## See also

- [`build_phyloseq()`](https://lmartinezgili.github.io/readyomics/reference/build_phyloseq.md)

- [`zCompositions::cmultRepl()`](https://rdrr.io/pkg/zCompositions/man/cmultRepl.html)

## Examples

``` r
if (requireNamespace("phyloseq", quietly = TRUE)) {
mock_X <- matrix(sample(0:1000, 25, replace = TRUE),
                 nrow = 5,
                 dimnames = list(paste0("sample", 1:5),
                 paste0("ASV", 1:5))
                 )

mock_sample_data <- data.frame(
  sample_id = paste0("sample", 1:5),
  load = c(1e5, 2e5, 1e4, 5e4, 1.5e5),
  condition = factor(rep(c("A", "B"), length.out = 5)),
  row.names = paste0("sample", 1:5)
  )

mock_taxa_table <- data.frame(
  Kingdom = rep("Bacteria", 5),
  Genus = paste0("Genus", 1:5),
  row.names = paste0("ASV", 1:5)
  )

result <- process_ngs(
  X = mock_X,
  sample_data = mock_sample_data,
  taxa_table = mock_taxa_table,
  normalise = "load",
  load_colname = "load",
  transform = "none",
  verbose = FALSE
  )
}
```
