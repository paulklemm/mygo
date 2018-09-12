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
    clusterProfiler::bitr(fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Mm.eg.db") %>%
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

#' Make GO-Term analysis and print out Report
#'
#' @export
#' @import magrittr ggplot2 rmarkdown
#' @param dat Dataframe containing variables `q_value` and `EntrezID` or `EnsemblID`
#' @param save_excel Save GO-terms to Excel files
#' @param significance_cutoff Specify the cutoff which entries are considered significant
#' @param output_path Path to HTML output
#' @return cummeRbund cuff object
createHTMLReport <- function(dat, output_path, save_excel = TRUE, significance_cutoff = 0.05) {
  # https://stackoverflow.com/questions/30377213/how-to-include-rmarkdown-file-in-r-package
  path_to_report <- system.file("rmd/Report.Rmd", package="mygo")
  # Render the document and put it into the output dir
  render(path_to_report, intermediates_dir = output_path, params = list(
    dat = dat,
    output_path = output_path,
    save_excel = save_excel,
    significance_cutoff = significance_cutoff
    ),
    output_dir = output_path,
    output_options=list(
      self_contained = TRUE
    )
  )
}

#' Create clusterProfiler EMA plot from go term
#'
#' @import magrittr clusterProfiler ggplot2
#' @param go_terms clusterProfiler enrichResult object
#' @param title ggtitle of the plot
ema_plot <- function(go_terms, title) {
  # There is a bug with very small GO term selections that we have to catch
  # HACK
  if (go_terms %>% nrow() > 2) {
    plot <- go_terms %>% enrichplot::emapplot() + ggplot2::ggtitle(title)  
    print(plot)
  }
}

#' Print a series of clusterProfiler plots
#'
#' @export
#' @import magrittr clusterProfiler ggplot2
#' @param go_terms clusterProfiler enrichResult object
plot_gos_normal <- function(go_terms) {
  if (go_terms %>% nrow() == 0) {
    print("No go terms to plot")
    return()
  }
  plot <- go_terms %>% barplot(showCategory=10) + ggplot2::ggtitle('Barplot of Top10 GO Terms')
  print(plot)
  go_terms %>% ema_plot('EMA Plot')
}

#' Perform a clusterPrfiler enrichGO analysis
#'
#' @export
#' @import magrittr clusterProfiler org.Mm.eg.db
#' @param ontoloty Can be either "BP", "CC", "MF"
#' @param enztezgenes List of entrezgenes to use for GO analysis
perform_enrichGO <- function(ontology, entrezgenes) {
  clusterProfiler::enrichGO(
      gene          = entrezgenes,
      OrgDb         = org.Mm.eg.db::org.Mm.eg.db,
      ont           = ontology,
      readable      = TRUE
    ) %>% 
    return()
}

#' Return a volcano plot with annotated top_n values
#'
#' @export
#' @import ggplot2 ggrepel
#' @param go_terms clusterProfiler enrichResult object
volcano_plot <- function(dat, label_top_n = 20) {
  dat %<>% mutate(Significant = pValue <= 0.05)
  plot <- dat %>%
    ggplot2::ggplot(aes(FoldChange, -log10(pValue))) +
    ggplot2::geom_point(aes(color = Significant)) +
    ggplot2::scale_color_manual(values = c("grey", "red")) +
    ggrepel::geom_text_repel(
        data = dat %>% filter(Significant == TRUE) %>% arrange(pValue) %>% head(label_top_n),
        aes(label = Symbol)
      )
  return(plot)
}

#' Get GO-terms for all ontologies
#'
#' @import magrittr dplyr
#' @export
#' @param dat Data frame containing columns `pValue` and `EntrezID`
#' @param significance_cutoff Cutoff to consider genes as significant
get_go_all_ontologies <- function(dat, significance_cutoff = 0.05) {
  significant_entrezgenes <- dat %>%
    dplyr::filter(q_value <= significance_cutoff) %>%
    dplyr::filter(!is.na(EntrezID)) %>%
    dplyr::mutate(EntrezID = as.numeric(EntrezID)) %>% .$EntrezID
  go_terms <- list()
  go_terms$BP <- significant_entrezgenes %>% perform_enrichGO(ontology = 'BP')
  go_terms$MF <- significant_entrezgenes %>% perform_enrichGO(ontology = 'MF')
  go_terms$CC <- significant_entrezgenes %>% perform_enrichGO(ontology = 'CC')
  go_terms %>% return()
}

#' Print a couple of standard plots for provided GO ontologies
#'
#' @export
#' @import purrr magrittr
#' @param ontologies List of clusterProfiler enrichResult objects with elements BP, MF and CC
plot_go_terms_all_ontologies <- function(ontologies) {
  ontologies %>%
  purrr::iwalk(function(.x, .y){
    # Print GO ontology
    .y %>% print()
    # Print GO ontology plots
    .x %>% plot_gos_normal()
  })
}

#' Export list of GO terms to Excel file
#'
#' @export
#' @import WriteXLS tibble magrittr
#' @param ontologies List of clusterProfiler enrichResult objects with elements BP, MF and CC
#' @param path Path to exported excel file
export_go_terms_to_excel <- function(ontologies, path) {
  BP <- ontologies$BP %>% as.tibble()
  MF <- ontologies$MF %>% as.tibble()
  CC <- ontologies$CC %>% as.tibble()
  c('BP', 'MF', 'CC') %>%
    WriteXLS::WriteXLS(ExcelFileName = path, AdjWidth = TRUE, AutoFilter = TRUE, BoldHeaderRow = TRUE, FreezeRow = 1, SheetNames = c('Biological Processes', 'Molecular Function', 'Cellular Component'))
}