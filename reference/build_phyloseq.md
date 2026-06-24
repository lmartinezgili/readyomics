# Build phyloseq objects for all taxonomy ranks

Constructs a list of phyloseq objects from a feature matrix (`X`),
sample data, taxonomy and (optionally) phylogenetic tree data.

## Usage

``` r
build_phyloseq(
  X,
  sample_data,
  taxa_table = NULL,
  phylo_tree = NULL,
  taxa_in_rows,
  verbose = TRUE
)
```

## Arguments

- X:

  A numeric matrix of NGS features (e.g., ASVs), with samples in rows
  and features in columns (recommended) or vice versa.

- sample_data:

  A `data.frame` containing sample data. Row names must match sample
  identifiers in `X`.

- taxa_table:

  (Optional) A taxonomy table with row names corresponding to feature
  names in `X`, and taxonomic ranks as columns.

- phylo_tree:

  (Optional) A phylogenetic tree.

- taxa_in_rows:

  Logical. If `TRUE`, `X` is assumed to have taxa as rows and samples as
  columns.

- verbose:

  Logical. If `TRUE`, diagnostic messages will be printed.

## Value

A named `list` of `phyloseq` objects and related output:

- asv:

  Phyloseq object with the raw feature counts (usually ASVs).

- \<tax_rank\>:

  Phyloseq objects of higher taxonomy ranks from `taxa_table`.

## Details

Phyloseq objects for higher taxonomic ranks are also generated when
`taxa_table` is provided. Higher rank taxa with labels matching
"unclass" or "unknown" are excluded after aggregation.

If very long strings are detected as feature IDs in `X` matrix or
`taxa_table`, (for example when actual DNA sequence is used as ID), it
will issue a warning, as this could significantly slow down computation
and increase memory usage.

## See also

[`phyloseq::phyloseq()`](https://rdrr.io/pkg/phyloseq/man/phyloseq.html)
for further details on phyloseq objects.

## Examples

``` r
if (requireNamespace("phyloseq", quietly = TRUE)) {
mock_X <- matrix(c(10, 0, 5, 3, 1, 7),
                 nrow = 2, byrow = TRUE,
                 dimnames = list(c("sample1", "sample2"),
                                 c("ASV1", "ASV2", "ASV3"))
                 )

mock_sample_data <- data.frame(sample_id = c("sample1", "sample2"),
                               group = c("A", "B"),
                               row.names = c("sample1", "sample2")
                               )

mock_taxa_table <- data.frame(Domain = c("Bacteria", "Bacteria", "Bacteria"),
                              Genus = c("GenusA", "GenusB", "Unknown"),
                              row.names = c("ASV1", "ASV2", "ASV3")
                              )

phyloseq_ready <- build_phyloseq(X = mock_X,
                                 sample_data = mock_sample_data,
                                 taxa_table = mock_taxa_table,
                                 taxa_in_rows = FALSE,
                                 verbose = FALSE)
}
```
