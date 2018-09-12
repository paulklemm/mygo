# mygo

Conduct GO-term analysis using [clusterProfiler](https://guangchuangyu.github.io/software/clusterProfiler/) and print report.

## Examples

```R
# Create data frame that fits the need of the analysis
dat <- readxl::read_xlsx('liver_chemogeentic.xlsx') %>%
  dplyr::rename(ensembl_gene_id = EnsemblID, q_value = pValue) %>%
  dplyr::select(ensembl_gene_id, q_value, Symbol)

dat %>% mygo::createHTMLReport(
  output_path = getwd()
)
```

## Restrictions & Data Preparations

Currently only `mus musculus` datasets are supported.

The data frame input needs to have the following dimensions:

- `ensembl_gene_id` or `EntrezID` (`character`)
- `q_value` (`numeric`)
- `Symbol` (`character`)

## Credits

- [clusterProfiler](https://guangchuangyu.github.io/software/clusterProfiler/) by Guangchuang Yu
