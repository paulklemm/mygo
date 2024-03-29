# mygo

---

<!-- TOC depthFrom:2 -->

- [Restrictions & Data Preparations](#restrictions--data-preparations)
- [Examples](#examples)
  - [Whole dataset](#whole-dataset)
  - [Selected genes](#selected-genes)
  - [Options for Rendering](#options-for-rendering)
- [Debug Run](#debug-run)
- [Installation](#installation)
- [Credits](#credits)
- [History](#history)

<!-- /TOC -->

---

Conduct GO-term analysis using [clusterProfiler](https://guangchuangyu.github.io/software/clusterProfiler/) and print report.

## Restrictions & Data Preparations

Currently only `mus musculus` datasets are supported.

The data frame input needs to have the following dimensions:

- `ensembl_gene_id` or `EntrezID` (`character`)
- `q_value` (`numeric`)
- `fc` (`numeric`)
- `Symbol` (`character`) (optional, but it's easier to identify genes that cannot be associated with an EntrezID)

Here is an example tibble.

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

## Examples

### Whole dataset

```R

library(magrittr)
# Create data frame that fits the need of the analysis
dat <- readr::read_tsv('test/geneexp_F_CPu.tsv') %>%
  dplyr::rename(fc = `log2(fold_change)`) %>%
  dplyr::mutate(Symbol = ensembl_gene_id) %>%
  dplyr::filter(status == "OK") %>%
  dplyr::select(ensembl_gene_id, q_value, fc)

# Here we do not pass a Symbol column. This won't result in an error.
dat %>% mygo::createHTMLReport(
  output_path = file.path(getwd(), 'result')
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
    do_gse = FALSE,
    use_background = TRUE
)

```

### Options for Rendering

- `dat` Input data frame
- `output_path` Output path of the analysis document and the excel file
- `save_excel`. Save GO-term result as Excel file
- `significance_cutoff`. Filter for significant GO terms
- `simplify_ontologies`. Run computational heavy GO term simplification.
- `do_gse`. Conduct a GSEA analysis. Deactivate if you do not pipe in a whole gene set.
- `use_background`. Use background of significant/not significant genes instead of all genes.

## Debug Run

Internal use only.

```r

# Local test run
my_dat <- readr::read_tsv('test/geneexp_F_CPu.tsv') %>%
  dplyr::rename(fc = `log2(fold_change)`) %>%
  dplyr::mutate(Symbol = ensembl_gene_id) %>%
  dplyr::filter(status == "OK") %>%
  dplyr::select(ensembl_gene_id, q_value, fc)

xaringan::infinite_moon_reader(
  moon = "inst/rmd/goterm_report.Rmd",
  cast_from=file.path(getwd(), "inst", "rmd"),
  params = list(
    dat = my_dat,
    output_path = ".",
    save_excel = FALSE,
    significance_cutoff = 0.05,
    simplify_ontologies = FALSE,
    do_gse = FALSE,
    use_background = TRUE,
    store_r_objects = FALSE,
    save_plots_as_pdf = FALSE
  )
)
```

Debug nfcore-rnaseq-pipeline run

```r

# Test for nfcore-rnaseq-pipeline output
nfcore_pipeline_path <- "/beegfs/scratch/bruening_scratch/pklemm/2020-11-anna-rnaseq/nfcore-rnaseq-pipeline"
dat <-
  glue::glue("{nfcore_pipeline_path}/results/DESeq2/npc/deseq_diff/deseq2_diff.csv") %>%
  readr::read_csv() %>%
  # Create data frame compatible with mygo
  dplyr::rename(
    ensembl_gene_id = row,
    q_value = padj,
    fc = log2FoldChange,
    Symbol = external_gene_name
  ) %>%
  dplyr::select(ensembl_gene_id, q_value, fc, Symbol)
xaringan::infinite_moon_reader(
  moon = "inst/rmd/goterm_report.Rmd",
  cast_from=file.path(getwd(), "inst", "rmd"),
  params = list(
    dat = dat,
    output_path = glue::glue("{nfcore_pipeline_path}/results/goterm-analyses/npc"),
    significance_cutoff = 0.1
  )
)

```


## Installation

```r
# install.packages('devtools')
devtools::install_github("paulklemm/mygo")
```

## Credits

- [clusterProfiler](https://guangchuangyu.github.io/software/clusterProfiler/) by Guangchuang Yu

## History

- *2023-07-07*
  - Add `log2FoldChange` filter to `run_goterms` and `get_go_all_ontologies`
  - Bump to 0.1.5
- *2022-09-30*
  - Add `run_kegg` option to `run_goterms`
  - Add `get_kegg_table` function
  - Bump to 0.1.4
- *2022-09-09*
  - Add `run_goterms` for usage in scripts
- *2022-02-10*
  - Add `get_go_all_ontologies_2` as a more flexible way of providing significant genes
  - Bump version to 0.1.1
- *2021-10-19*
  - Add `bind_goterm_table`, `emap_plot_facet_category`, `overlap_percentage_plot_facet_category`
  - Remove all imports except `magrittr`
  - Bump version to 0.1.0
- *2021-06-29*
  - Add `set_readable` of DOSE to fix ridgeplot
  - Add `simplify` option for GSEA
  - Bump to `0.0.15`
- *2021-04-12*
  - Fix sorting of `overlap_percentage_plot`
  - Bump version to `0.0.13`
- *2020-11-17*
  - Fix `emap_plot`
  - Bump version to `0.0.12`
- *2020-10-26*
  - Add `print_kegg`
  - Bump version to `0.0.10`
- *2020-09-09*
  - Export `print_goterm_as_datatable` function
  - Bump version to `0.0.9`
- *2020-09-07*
  - Add data dictionary
- *2020-05-12*
  - Knit the document in a temporary folder to avoid files being cleaned up in the render function
- *2020-04-22*
  - Have complete tables in the report file so users can download their tables from there
- *2020-03-25*
  - Remove the overlap percentage in favor or plotting the top-regulated terms
- *2020-02-18*
  - Remove `radix`, add `ggthemes` to `DESCRIPTION`
- *2020-02-07*
  - Added `save_plots_as_pdf` option
- *2020-02-05-06*
  - Completely overhaul mygo output
- *2020-02-04*
  - Added `store_r_objects` option
  - Added `simplify_cutoff` option. See [https://guangchuangyu.github.io/2015/10/use-simplify-to-remove-redundancy-of-enriched-go-terms/](https://guangchuangyu.github.io/2015/10/use-simplify-to-remove-redundancy-of-enriched-go-terms/)
