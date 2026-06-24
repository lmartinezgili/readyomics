# Add taxonomic information to `dana` object

Appends features taxonomy to the `dana` object tables.

## Usage

``` r
add_taxa(
  dana_obj,
  taxa_table,
  taxa_rank = c("asv", "substrain", "strain", "species", "genus", "family", "order",
    "class", "phylum", "domain")
)
```

## Arguments

- dana_obj:

  A `dana` object returned by
  [`dana()`](https://lmartinezgili.github.io/readyomics/reference/dana.md).

- taxa_table:

  A taxonomy table `data.frame` with taxonomy ranks in columns and row
  names corresponding to `feat_id`s in `dana` object.

- taxa_rank:

  A character string specifying the taxonomy level of input features.
  Accepts one of: `"asv"`, `"substrain"`, `"strain"`, `"species"`,
  `"genus"`, `"family"`, `"order"`, `"class"`, `"phylum"`, or
  `"domain"`.

## Value

A modified version of `dana_obj`, with taxonomy information added to
relevant tables.

## Details

- If `taxa_rank = "asv"`, a `taxon_name` is constructed by pasting the
  ASV ID to the `species` (if available) or `genus` name.

- For other ranks, `taxon_name` is taken directly from the corresponding
  column in `taxa_table`.

- All higher-level taxonomy ranks available in `taxa_table` are also
  appended.

## See also

[`dana()`](https://lmartinezgili.github.io/readyomics/reference/dana.md)
for fitting differential analysis models on omics datasets.

## Examples

``` r
set.seed(123)
mock_X <- matrix(rnorm(20 * 5), nrow = 20)
colnames(mock_X) <- paste0("feat_", seq_len(5))
rownames(mock_X) <- paste0("sample_", seq_len(20))

mock_taxa <- data.frame(
  Domain = rep("Bacteria", 5),
  Phylum = c("Firmicutes", "Bacteroidota", "Proteobacteria",
             "Actinobacteriota", "Firmicutes"),
  Class = c("Bacilli", "Bacteroidia", "Gammaproteobacteria",
            "Actinobacteria", "Clostridia"),
  Order = c("Lactobacillales", "Bacteroidales", "Enterobacterales",
            "Bifidobacteriales", "Clostridiales"),
  Family = c("Lactobacillaceae", "Bacteroidaceae", "Enterobacteriaceae",
             "Bifidobacteriaceae", "Clostridiaceae"),
  Genus = c("Lactobacillus", "Bacteroides", "Escherichia",
            "Bifidobacterium", "Clostridium"),
  Species = c("acidophilus", "fragilis", "coli", "longum", "butyricum"),
  row.names = paste0("feat_", seq_len(5)),
  stringsAsFactors = FALSE
)

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
  padj = runif(10),
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

# Add taxonomy
dana_obj <- dana_obj |>
  add_taxa(mock_taxa, taxa_rank = "genus")
```
