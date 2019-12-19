# mygo

---

<!-- TOC depthFrom:2 -->

- [Examples](#examples)
  - [Whole dataset](#whole-dataset)
  - [Selected genes](#selected-genes)
  - [Options for Rendering](#options-for-rendering)
- [Installation](#installation)
- [Restrictions & Data Preparations](#restrictions--data-preparations)
- [Credits](#credits)

<!-- /TOC -->

---

Conduct GO-term analysis using [clusterProfiler](https://guangchuangyu.github.io/software/clusterProfiler/) and print report.

## Examples

MyGO expects a data frame containing `ensembl_gene_id`/`EntrezID`, `q_value` and `fc`.
If `EntrezID` is not included, it will convert `ensembl_gene_id` to `EntrezID`

```R
# A tibble: 18,777 x 3
   ensembl_gene_id    q_value       fc
   <chr>                <dbl>    <dbl>
 1 ENSMUSG00000103922  0.998   0.00361
 2 ENSMUSG00000025903  0.443   0.141  
 3 ENSMUSG00000104217  0.443   0.141  
 4 ENSMUSG00000033813  0.443   0.141  
 5 ENSMUSG00000033793  0.718   0.0864 
 6 ENSMUSG00000025905  0.0553  0.334  
 7 ENSMUSG00000025907  0.747   0.0919 
 8 ENSMUSG00000087247  0.282  -3.24   
 9 ENSMUSG00000033740  0.431  -0.221  
10 ENSMUSG00000102135  0.672  -0.547  
# … with 18,767 more rows
```

### Whole dataset

```R
library(magrittr)
# Create data frame that fits the need of the analysis
dat <- readr::read_tsv('test/geneexp_F_CPu.tsv') %>%
  dplyr::rename(fc = `log2(fold_change)`) %>%
  dplyr::mutate(Symbol = ensembl_gene_id) %>%
  dplyr::filter(status == "OK") %>%
  dplyr::select(ensembl_gene_id, q_value, fc)
# Standard Mode
dat %>% mygo::createHTMLReport(
  output_path = file.path(getwd(), 'result')
)
# Debug Mode
dat %>% mygo::createHTMLReport(
  output_path = file.path(getwd(), 'result'),
  dev = TRUE
)
```

### Selected genes

Sometimes you want to do a GO term analysis only for a small number of genes. In this case, we need to make sure to get the proper gene names and deactivate the GSEA analysis.

Here is an example call.

```R
my_genes <- readxl::read_xlsx("data/dat.xlsx") %>%
  dplyr::mutate(
    gene_name = `Gene names`,
    fc = `-Log t-test p value`) %>%
  dplyr::select(
    gene_name,
    fc
  ) %>%
  tidyr::separate_rows(gene_name, sep = ";") %>%
  dplyr::distinct(gene_name, .keep_all = TRUE) %>%
  dplyr::mutate(
    gene_name_fixed = rmyknife::get_gene_name_from_synonym(gene_name)
  ) %>%
  rmyknife::attach_ensembl_gene_id_from_name(
    gene_name_var = "gene_name_fixed",
    ensembl_version = 96
  ) %>%
  # We do not need the q-values, so set them to 0
  dplyr::mutate(q_value = 0)

my_genes %>%
  mygo::createHTMLReport(
    output_path = file.path(getwd(), 'result', 'c5'),
    save_excel = TRUE,
    do_gse = FALSE
    use_background = TRUE
)
```


### Options for Rendering

- `dat` Input data frame
- `output_path` Output path of the analysis document and the excel file
- `save_excel`. Save GO-term result as Excel file
- `significance_cutoff`. Filter for significant GO terms
- `dev`. Developer mode. Do not use this.
- `simplify_ontologies`. Run computational heavy GO term simplification.
- `do_gse`. Conduct a GSEA analysis. Deactivate if you do not pipe in a whole gene set.

## Installation

```r
# install.packages('devtools')
devtools::install_github("paulklemm/mygo")
```

## Restrictions & Data Preparations

Currently only `mus musculus` datasets are supported.

The data frame input needs to have the following dimensions:

- `ensembl_gene_id` or `EntrezID` (`character`)
- `q_value` (`numeric`)
- `Symbol` (`character`)

## Credits

- [clusterProfiler](https://guangchuangyu.github.io/software/clusterProfiler/) by Guangchuang Yu
