# mygo

<!-- TOC depthFrom:2 -->

- [Examples](#examples)
- [Installation](#installation)
- [Restrictions & Data Preparations](#restrictions--data-preparations)
- [Credits](#credits)

<!-- /TOC -->

Conduct GO-term analysis using [clusterProfiler](https://guangchuangyu.github.io/software/clusterProfiler/) and print report.

## Examples

```R
library(tidyverse)
# Create data frame that fits the need of the analysis
dat <- readr::read_tsv('test/geneexp_F_CPu.tsv') %>%
  dplyr::rename(fc = `log2(fold_change)`) %>%
  dplyr::mutate(Symbol = ensembl_gene_id) %>%
  dplyr::filter(status == "OK") %>%
  dplyr::select(ensembl_gene_id, q_value, fc)

dat %>% mygo::createHTMLReport(
  output_path = getwd()
)
```

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
