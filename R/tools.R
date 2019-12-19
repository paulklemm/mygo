#' Gets Entrez IDs and add it to dat
#'
#' @export
#' @import magrittr clusterProfiler org.Mm.eg.db dplyr tibble
#' @param dat Data frame with ensembl identifier
#' @param ensembl_id_name Name of column containing the ensembl identifier
#' @param keep_only_rows_with_entrez Only keep rows for which entrez IDs could be found
#' @param drop_duplicates Often there is a n:1 relationship between Entrez-IDs and Ensembl-IDs. If this value is true, only keep the first hit
ensembl_to_entrez <- function(dat, ensembl_id_name, keep_only_rows_with_entrez = TRUE, drop_duplicates = TRUE) {
  # Get the Entrez IDs
  ens_to_ent <- dat[ensembl_id_name][[1]] %>%
    clusterProfiler::bitr(fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db") %>%
    dplyr::rename(EntrezID = ENTREZID)
  # Drop duplicates if required
  if (drop_duplicates) {
    ens_to_ent <- ens_to_ent %>% dplyr::filter(!duplicated(ENSEMBL))
  }
  if (keep_only_rows_with_entrez) {
    dat %>%
      dplyr::right_join(ens_to_ent, by = setNames("ENSEMBL", ensembl_id_name)) %>%
      tibble::as.tibble() %>%
      return()
  } else {
    dat %>%
      dplyr::full_join(ens_to_ent, by = setNames("ENSEMBL", ensembl_id_name)) %>%
      tibble::as.tibble() %>%
      return()
  }
}

#' Add Gene Symbol from EntrezID
#'
#' @export
#' @import clusterProfiler dplyr magrittr
#' @param dat data frame with EnzrezID
attach_gene_symbol_from_entrez <- function(dat) {
  if (!('EntrezID' %in% colnames(dat))) {
    stop('EntrezID not available in data frame')
  }

  dat["EntrezID"][[1]] %>% clusterProfiler::bitr(fromType = "ENTREZID", toType = c("SYMBOL"), OrgDb = org.Mm.eg.db) %>%
    dplyr::rename(EntrezID = ENTREZID, Symbol = SYMBOL) %>%
    dplyr::right_join(dat) %>%
    as.tibble() %>%
    return()
}

#' Make GO-Term analysis and print out Report
#'
#' @export
#' @import magrittr ggplot2 rmarkdown
#' @param dat Dataframe containing variables `q_value` and `EntrezID` or `EnsemblID`
#' @param save_excel Save GO-terms to Excel files
#' @param significance_cutoff Specify the cutoff which entries are considered significant
#' @param output_path Path to HTML output
#' @param simplify_ontologies Run time-consuming ontology simplification
#' @param do_gse Conduct GSE analysis
#' @param use_background Use background genes for the analysis
#' @return cummeRbund cuff object
createHTMLReport <- function(dat, output_path, save_excel = TRUE, significance_cutoff = 0.05, do_gse = TRUE, simplify_ontologies = TRUE, use_background = TRUE) {
  # https://stackoverflow.com/questions/30377213/how-to-include-rmarkdown-file-in-r-package
  path_to_report <- system.file("rmd/Report.Rmd", package = "mygo")
  # Render the document and put it into the output dir
  render(
    path_to_report,
    params = list(
      dat = dat,
      output_path = output_path,
      save_excel = save_excel,
      significance_cutoff = significance_cutoff,
      simplify_ontologies = simplify_ontologies,
      do_gse = do_gse,
      use_background = use_background
    ),
    # RMarkdown options
    output_dir = output_path,
    output_options = list(
      self_contained = TRUE
    ),
    intermediates_dir = output_path,
    knit_root_dir = output_path,
    clean = TRUE
  )
}

#' Create clusterProfiler EMA plot from go term
#'
#' @import magrittr clusterProfiler ggplot2 enrichplot
#' @param go_terms clusterProfiler enrichResult object
#' @param title ggtitle of the plot
emap_plot <- function(go_terms, title) {
  # There is a bug with very small GO term selections that we have to catch
  # HACK
  if (go_terms %>% nrow() > 5) {
    plot <- go_terms %>% enrichplot::emapplot() + ggplot2::ggtitle(title)
    return(plot)
  }
}

#' Print a series of clusterProfiler plots
#'
#' @export
#' @import magrittr clusterProfiler ggplot2 enrichplot plotly
#' @param go_terms clusterProfiler enrichResult object
#' @param fc_symbol double vector containing the fold changes named by symbol
plot_terms_go <- function(go_terms, fc_symbol) {
  if (go_terms %>% nrow() == 0) {
    print("No go terms to plot")
    return()
  }
  # go_terms %>% clusterProfiler::cnetplot(foldChange = fc_symbol, circular = TRUE, colorEdge = TRUE) %>% print()
  # plot <- go_terms %>% barplot(showCategory=10) + ggplot2::ggtitle('Barplot of Top10 GO Terms')
  # plotly::ggplotly(plot)
  plot <- go_terms %>% barplot(showCategory = 10) + ggplot2::ggtitle('Barplot of Top10 GO Terms')
  plot %>% print()
  # We had a problem getting the GOplot for small numbers of GO-terms
  tryCatch(
    go_terms %>% enrichplot::goplot() %>% print(),
    error = function(error_message) { return("GOPlot not available") }
  )
  go_terms %>% emap_plot('Enrich Map Plot Plot') %>% print()
}

#' Plot GSE terms
#'
#' @export
#' @import clusterProfiler ggplot2
#' @param gse_terms clusterProfiler gse terms object
#' @param fc_symbol double vector containing the fold changes named by symbol
plot_terms_gse <- function(gse_terms, fc_symbol) {
  if (gse_terms %>% nrow() == 0) {
    print("No go terms to plot")
    return()
  }
  gse_terms %>%
    enrichplot::emapplot() %>%
    print()
  gse_terms %>%
    enrichplot::heatplot(foldChange = fc_symbol) %>%
    print()
  #enrichplot::heatplot(foldChange = fc_symbol) %>%
  #plotly::ggplotly()
  gse_terms %>%
    enrichplot::gseaplot(geneSetID = 1) %>%
    print()
  gse_terms %>%
    enrichplot::cnetplot(foldChange = fc_symbol, circular = TRUE, colorEdge = TRUE) %>%
    print()
  gse_terms %>%
    enrichplot::ridgeplot() %>%
    print()
}

#' Simplify GO terms
#'
#' @export
#' @import clusterProfiler magrittr purrr
#' @param ontologies List of clusterProfiler enrichResult objects with elements BP, MF and CC
simplify_ontologies <- function(ontologies) {
  ontologies %>%
  purrr::map(function(.x) {
    # Print GO ontology plots
    .x %>% clusterProfiler::simplify()
  }) %>%
    return()
}

#' Perform a clusterPrfiler enrichGO analysis
#'
#' @export
#' @import magrittr clusterProfiler org.Mm.eg.db
#' @param ontology Can be either "BP", "CC", "MF"
#' @param entrezgenes List of entrezgenes to use for GO analysis
#' @param entrez_background_genes List of background genes 
#' @param use_background use specified background genes
perform_enrichGO <- function(ontology, entrezgenes, background_genes, use_background) {
  if (use_background) {
    clusterProfiler::enrichGO(
      gene = entrezgenes,
      OrgDb = org.Mm.eg.db::org.Mm.eg.db,
      universe = background_genes,
      ont = ontology,
      readable = TRUE
    ) %>%
    return()
  } else {
    clusterProfiler::enrichGO(
      gene = entrezgenes,
      OrgDb = org.Mm.eg.db::org.Mm.eg.db,
      ont = ontology,
      readable = TRUE
    ) %>%
    return()
  }
}

#' Return a volcano plot with annotated top_n values
#'
#' @export
#' @import ggplot2 ggrepel
#' @param go_terms clusterProfiler enrichResult object
volcano_plot <- function(dat, label_top_n = 20) {
  dat %<>% mutate(Significant = q_value <= 0.05)
  plot <- dat %>%
    ggplot2::ggplot(aes(fc, - log10(q_value), key = Symbol)) +
    ggplot2::geom_point(aes(color = Significant)) +
    ggplot2::scale_color_manual(values = c("grey", "red"))
  return(plotly::ggplotly(plot))
}

#' Get GO-terms for all ontologies
#'
#' @import magrittr dplyr
#' @export
#' @param dat Data frame containing columns `pValue` and `EntrezID`
#' @param significance_cutoff Cutoff to consider genes as significant
#' @param use_background Use gene list as background instead of using all genes as background
get_go_all_ontologies <- function(dat, significance_cutoff = 0.05, use_background = TRUE) {
  # Prepare data frame
  valid_dat <- dat %>%
    dplyr::filter(!is.na(EntrezID)) %>%
    #dplyr::mutate(EntrezID = as.numeric(EntrezID))
    dplyr::mutate(EntrezID = as.character(EntrezID))

  # Get all genes in the data set as background universe
  background_entrezgenes <- valid_dat %>%
    .$EntrezID
  # Get significant genes
  significant_entrezgenes <- valid_dat %>%
    dplyr::filter(q_value <= significance_cutoff) %>%
    .$EntrezID
  # Perform GO-term analysis for each term
  go_terms <- list()
  go_terms$Biological_Process <- significant_entrezgenes %>% perform_enrichGO(
    ontology = "BP",
    background_genes = background_entrezgenes,
    use_background = use_background
  )
  go_terms$Molecular_Function <- significant_entrezgenes %>% perform_enrichGO(
    ontology = "MF",
    background_genes = background_entrezgenes,
    use_background = use_background
  )
  go_terms$Cellular_Components <- significant_entrezgenes %>% perform_enrichGO(
    ontology = "CC",
    background_genes = background_entrezgenes,
    use_background = use_background
  )
  go_terms %>% return()
}

#' Print a couple of standard plots for provided GO ontologies
#'
#' @export
#' @import purrr magrittr
#' @param ontologies List of clusterProfiler enrichResult objects with elements BP, MF and CC
#' @param fc_symbol double vector containing the fold changes named by symbol
plot_all_ontologies <- function(ontologies, fc_symbol) {
  ontologies %>%
  purrr::iwalk(function(.x, .y) {
    # Print GO ontology
    .y %>% print()
    # Print GO ontology plots
    if (.x %>% .hasSlot("setType")) {
      .x %>% plot_terms_gse(fc_symbol)
    } else {
      .x %>% plot_terms_go(fc_symbol)
    }
  })
}

#' Export list of GO terms to Excel file
#' TODO: Redo this!
#'
#' @export
#' @import WriteXLS tibble magrittr
#' @param go_ontologies List of clusterProfiler enrichResult objects with elements BP, MF and CC
#' @param go_ontologies_simple List of clusterProfiler enrichResult objects with elements BP, MF and CC after being simplified
#' @param gse_ontologies GSE enrichResult
#' @param kegg_ontologies KEGG enrichResult
#' @param path Path to exported excel file
export_go_terms_to_excel <- function(go_ontologies, go_ontologies_simple, gse_ontologies, kegg_ontologies, path) {
  BP_go <- go_ontologies$Biological_Process %>% as.tibble()
  MF_go <- go_ontologies$Molecular_Function %>% as.tibble()
  CC_go <- go_ontologies$Cellular_Components %>% as.tibble()
  BP_go_simple <- go_ontologies_simple$Biological_Process %>% as.tibble()
  MF_go_simple <- go_ontologies_simple$Molecular_Function %>% as.tibble()
  CC_go_simple <- go_ontologies_simple$Cellular_Components %>% as.tibble()
  BP_gse <- gse_ontologies$Biological_Process %>% as.tibble()
  MF_gse <- gse_ontologies$Molecular_Function %>% as.tibble()
  CC_gse <- gse_ontologies$Cellular_Components %>% as.tibble()
  kegg_gse <- kegg_ontologies$kegg %>% as.tibble()

  c('BP_go', 'MF_go', 'CC_go', 'BP_go_simple', 'MF_go_simple', 'CC_go_simple', 'BP_gse', 'MF_gse', 'CC_gse', 'kegg_gse') %>%
    WriteXLS::WriteXLS(ExcelFileName = path,
      AdjWidth = TRUE,
      AutoFilter = TRUE,
      BoldHeaderRow = TRUE,
      FreezeRow = 1,
      SheetNames = c('GO Biological Processes', 'GO Molecular Function', 'GO Cellular Component', 'GO Simple Biological Processes', 'GO Simple Molecular Function', 'GO Simple Cellular Component', 'GSE Biological Processes', 'GSE Molecular Function', 'GSE Cellular Component', 'GSE KEGG')
    )
}

#' Get named fc vector. Prepare input data as required by GSE
#'
#' @export
#' @import magrittr clusterProfiler
#' @param dat Dataframe with columns `EntrezID` and `fc`
get_named_fc_vector <- function(dat) {
  # https://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html#go-gene-set-enrichment-analysis
  fc <- dat %>% .$fc %>% as.double() %>% `names<-`(dat$EntrezID)
  fc <- fc[!is.infinite(fc)]
  fc <- fc %>% sort(decreasing = TRUE)
  fc %>% return()
}

#' Get GSE-terms for all ontologies
#'
#' @import magrittr dplyr
#' @export
#' @param dat Data frame containing columns `pValue` and `EntrezID`
get_gse_all_ontologies <- function(dat) {
  # Prepare input data as required by GSE
  fc <- dat %>% get_named_fc_vector()
  gse_terms <- list()
  # Barplot of the fold change
  # fc %>% barplot()
  gse_terms$Biological_Process <- perform_gseGO(ontology = 'BP', fc = fc)
  gse_terms$Molecular_Function <- perform_gseGO(ontology = 'MF', fc = fc)
  gse_terms$Cellular_Components <- perform_gseGO(ontology = 'CC', fc = fc)
  gse_terms %>% return()
}

#' Get KEGG GSE terms
#'
#' @import magrittr dplyr
#' @export
#' @param dat Data frame containing columns `pValue` and `EntrezID`
get_kegg <- function(dat) {
  # Prepare input data as required by GSE
  fc <- dat %>% get_named_fc_vector()
  kegg_terms <- list()
  kegg_terms$kegg <- perform_gseKEGG(fc = fc)
  kegg_terms %>% return()
}

#' Perform a clusterProfiler gseGO analysis
#'
#' @export
#' @import magrittr clusterProfiler org.Mm.eg.db DOSE
#' @param ontoloty Can be either "BP", "CC", "MF"
#' @param fc Named vector of foldchanges (name denotes Entrez ID)
perform_gseGO <- function(ontology, fc) {
  clusterProfiler::gseGO(geneList = fc,
    OrgDb = org.Mm.eg.db,
    ont = ontology,
    nPerm = 1000,
    minGSSize = 100,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    verbose = FALSE
  ) %>%
    DOSE::setReadable(OrgDb = org.Mm.eg.db) %>%
    return()
}

#' Perform a clusterProfiler gseKEGG analysis
#'
#' @export
#' @import magrittr clusterProfiler DOSE
#' @param fc Named vector of foldchanges (name denotes Entrez ID)
perform_gseKEGG <- function(fc) {
  fc %>% clusterProfiler::gseKEGG(nPerm = 10000, organism = 'mmu') %>%
    return()
}
