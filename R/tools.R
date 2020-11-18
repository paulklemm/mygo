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
#' @param store_r_objects Store clusterProfiler R objects as compressed RDS files to output_path
#' @param simplify_cutoff See clusterProfiler::simplify
#' @param save_plots_as_pdf Export all plots as PDF to output_path/plots
#' @return cummeRbund cuff object
createHTMLReport <- function(
  dat,
  output_path,
  save_excel = TRUE,
  significance_cutoff = 0.05,
  do_gse = TRUE,
  simplify_ontologies = TRUE,
  use_background = TRUE,
  store_r_objects = TRUE,
  simplify_cutoff = 0.7,
  save_plots_as_pdf = TRUE
  ) {
  # https://stackoverflow.com/questions/30377213/how-to-include-rmarkdown-file-in-r-package
  path_to_report <- system.file("rmd/goterm_report.Rmd", package = "mygo")
  # Create temporary knitting dir
  temp_knitting_path <- file.path(output_path, "knitting")
  dir.create(path = temp_knitting_path, recursive = TRUE)
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
      use_background = use_background,
      store_r_objects = store_r_objects,
      simplify_cutoff = simplify_cutoff,
      save_plots_as_pdf = save_plots_as_pdf
    ),
    # RMarkdown options
    output_dir = output_path,
    output_file = "goterm_report.html",
    output_options = list(
      self_contained = TRUE
    ),
    intermediates_dir = temp_knitting_path,
    knit_root_dir = temp_knitting_path,
    clean = TRUE
  )
  # Now remove the intermediate knitting directory
  file.remove(temp_knitting_path)
}

#' Create clusterProfiler EMA plot from go term
#'
#' @import magrittr clusterProfiler ggplot2 enrichplot
#' @export
#' @param go_terms clusterProfiler enrichResult object
#' @param title ggtitle of the plot
#' @param n Number of GO-term to plot
emap_plot <- function(go_terms, title, n = 30) {
  # Catch null emap plot
  if (is.null(go_terms)) {
    warning("Cannot print emap plot, no GO-terms")
    return()
  }
  # There is a bug with very small GO term selections that we have to catch
  # HACK
  if (go_terms %>% nrow() <= 5) {
    return()
  }
  
  # Get myplot using tryCatch
  myplot <- tryCatch(
    expr = {
      myplot <- go_terms %>%
        # Comply with latest implementation https://github.com/YuLab-SMU/clusterProfiler/issues/299
        enrichplot::pairwise_termsim() %>%
        # Run emapplot
        enrichplot::emapplot(showCategory = n) +
        ggplot2::ggtitle(title)
      return(myplot)
    },
    error = function(e){
      message("Can't print emap. Possibly too few differentially expressed terms. Error:")
      message(e)
      return(NULL)
    }
  )
  return(myplot)
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
  # There is a bug with very small GO term selections that we have to catch
  # HACK
  if (gse_terms %>% nrow() > 5) {
    gse_terms %>%
      # Comply with latest implementation https://github.com/YuLab-SMU/clusterProfiler/issues/299
      enrichplot::pairwise_termsim() %>%
      enrichplot::emapplot() %>%
      print()
  }
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
#' @param cutoff See clusterProfiler::simplify
simplify_ontologies <- function(ontologies, cutoff) {
  ontologies %>%
  purrr::map(function(.x) {
    # Print GO ontology plots
    .x %>% clusterProfiler::simplify(cutoff = cutoff)
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


#' Add count for each GO-term and calculate overlap ratio of significant genes
#' @export
#' @import dplyr magrittr tibble
#' @param dat clusterProfiler GO-term result DOSE
attach_goterm_genecount <- function(
  dat
) {
  dat_tibble <- dat %>%
    tibble::as_tibble()
  # Check if we have a populated dataframe
  if (dat_tibble %>% nrow() > 0) {
    dat_tibble %>%
      dplyr::rowwise() %>%
      tidyr::separate(GeneRatio, remove = FALSE, into = c("temp", "TotalCount")) %>%
      # Remove temporary variable
      dplyr::select(-temp) %>%
      tidyr::separate(BgRatio, remove = FALSE, into = c("GOTermGeneCount", "BackgroundCount")) %>%
      # Convert variables to numeric
      dplyr::mutate_at(c("TotalCount", "GOTermGeneCount", "BackgroundCount"), as.numeric) %>%
      dplyr::mutate(
        # Do not calculate GO-term size from all genes, but rather from the background
        # GOTermGeneCount = rmyknife::get_genes_of_goterm_godb(ID, species = species, verbose = FALSE) %>% length(),
        # This is also called a rich_factor
        Percent_Significant = (Count * 100) / GOTermGeneCount
        # See https://yulab-smu.github.io/clusterProfiler-book/chapter13.html
        # FoldEnrichment = (Count / TotalCount) / (GOTermGeneCount / BackgroundCount)
      ) %>%
      return()
  } else {
    warning("Cannot attach metrics to empty data frame")
    return(dat_tibble)
  }
}

#' Export list of GO terms to Excel file
#'
#' @export
#' @import WriteXLS tibble magrittr
#' @param go_ontologies List of clusterProfiler enrichResult objects with elements BP, MF and CC
#' @param go_ontologies_simple List of clusterProfiler enrichResult objects with elements BP, MF and CC after being simplified
#' @param gse_ontologies GSE enrichResult
#' @param kegg_ontologies KEGG enrichResult
#' @param path Path to exported excel file
#' @param species Either "MUS" or "HUM"
export_go_terms_to_excel <- function(
  go_ontologies,
  go_ontologies_simple,
  gse_ontologies,
  kegg_ontologies,
  path,
  species = "MUS"
) {
  BP_go <- go_ontologies$Biological_Process %>% attach_goterm_genecount()
  MF_go <- go_ontologies$Molecular_Function %>% attach_goterm_genecount()
  CC_go <- go_ontologies$Cellular_Components %>% attach_goterm_genecount()
  BP_go_simple <- go_ontologies_simple$Biological_Process %>% attach_goterm_genecount()
  MF_go_simple <- go_ontologies_simple$Molecular_Function %>% attach_goterm_genecount()
  CC_go_simple <- go_ontologies_simple$Cellular_Components %>% attach_goterm_genecount()
  # GSE objects have no GeneRatio object, so we cannot add the goterm genecount
  BP_gse <- gse_ontologies$Biological_Process %>% tibble::as_tibble()
  MF_gse <- gse_ontologies$Molecular_Function %>% tibble::as_tibble()
  CC_gse <- gse_ontologies$Cellular_Components %>% tibble::as_tibble()
  kegg_gse <- kegg_ontologies$kegg %>% tibble::as_tibble()

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
#' @param p_cutoff P-value cutoff for GSEA
get_gse_all_ontologies <- function(dat, p_cutoff = 0.05) {
  # Prepare input data as required by GSE
  fc <- dat %>% get_named_fc_vector()
  gse_terms <- list()
  # Barplot of the fold change
  # fc %>% barplot()
  gse_terms$Biological_Process <- perform_gseGO(ontology = 'BP', fc = fc, p_cutoff = p_cutoff)
  gse_terms$Molecular_Function <- perform_gseGO(ontology = 'MF', fc = fc, p_cutoff = p_cutoff)
  gse_terms$Cellular_Components <- perform_gseGO(ontology = 'CC', fc = fc, p_cutoff = p_cutoff)
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
#' @param p_cutoff P-value cutoff for GSEA
perform_gseGO <- function(ontology, fc, p_cutoff = 0.05) {
  clusterProfiler::gseGO(
    geneList = fc,
    OrgDb = org.Mm.eg.db,
    ont = ontology,
    pvalueCutoff = p_cutoff,
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
  fc %>% clusterProfiler::gseKEGG(organism = "mmu") %>%
    return()
}

#' Plot total gene count against percentage of significant genes
#' @export
#' @import ggrepel ggplot2 magrittr dplyr ggthemes
#' @param dat clusterProfiler GO-term table with attached GOTermGeneCount and Percent_Significant column
#' @param n Number of GO-terms to label
#' @examples
#'   dat_goterms$Biological_Process %>%
#'     mygo::attach_goterm_genecount() %>%
#'     mygo::overlap_scatterplot
overlap_scatterplot <- function(
  dat,
  n = 15
) {
  if (dat %>% nrow() == 0) {
    warning("Cannot create overlap scatterplot for empty data frame")
    return()
  }
  plot <- dat %>%
    ggplot2::ggplot(
      mapping = ggplot2::aes(
        x = GOTermGeneCount,
        y = Percent_Significant,
        color = p.adjust
      )
    ) +
    ggplot2::geom_point(alpha = 0.6) +
    ggplot2::scale_color_continuous(high = "#132B43", low = "#56B1F7") +
    ggrepel::geom_label_repel(
      data = . %>% dplyr::arrange(p.adjust) %>% head(n),
      mapping = ggplot2::aes(label = Description),
      color = "black",
      min.segment.length = 0
    ) +
    ggplot2::theme_minimal() +
    ggplot2::xlab("Total GO-term gene count") + 
    ggplot2::ylab("Percentage of significant genes")
  return(plot)
}

#' Plot top n GO-terms w.r.t percentage of significant genes
#' Plot adapted from https://yulab-smu.github.io/clusterProfiler-book/chapter13.html
#' @param dat clusterProfiler GO-term table with attached GOTermGeneCount and Percent_Significant column
#' @param n Number of GO-terms to label
#' @order_by Order by overlap "percentage" or "significance"
#' @export
#' @import ggplot2 forcats magrittr
#' @examples
#'   dat_goterms$Biological_Process %>%
#'     mygo::attach_goterm_genecount() %>%
#'     mygo::overlap_percentage_plot()
overlap_percentage_plot <- function(dat, n = 25, order_by = "percentage") {
  if (dat %>% nrow() == 0) {
    warning("Cannot create overlap percentage plot for empty data frame")
    return()
  }
  # Order dat in the required way
  plot_title <- ""
  if (order_by == "percentage") {
    dat %<>% dplyr::arrange(dplyr::desc(Percent_Significant))
    plot_title <- paste0("Top ", n, " terms with highest overlap")
  } else if (order_by == "significance") {
    dat %<>% dplyr::arrange(p.adjust)
    plot_title <- paste0("Top ", n, " sig. terms, ordered by overlap")
  } else {
    warning(paste0("I don't know what to do with order_by parameter value ", order_by))
    return()
  }
  plot <- dat %>%
    dplyr::rename(`Significant Gene Count` = Count) %>%
    head(n) %>%
    ggplot2::ggplot(
      mapping = ggplot2::aes(
        Percent_Significant,
        forcats::fct_reorder(Description, Percent_Significant)
      )
    ) +
      ggplot2::geom_segment(ggplot2::aes(xend = 0, yend = Description)) +
      ggplot2::geom_point(ggplot2::aes(color = p.adjust, size = `Significant Gene Count`)) +
      # ggplot2::scale_color_viridis_c(guide = ggplot2::guide_colorbar(reverse = TRUE)) +
      ggplot2::scale_color_continuous(high = "#132B43", low = "#56B1F7") +
      ggplot2::scale_size_continuous(range = c(1, 5)) +
      ggplot2::theme_minimal() +
      ggplot2::xlab("Percentage of Significant Genes") +
      ggplot2::ylab(NULL) +
      ggplot2::ggtitle(plot_title)
  return(plot)
}

#' Helper function, takes clusterProfiler GO-term result and prints it as DT widget
#' @export
#' @import DT dplyr tibble magrittr
#' @param dat clusterProfiler GO-term result
print_goterm_as_datatable <- function(dat) {
  # Check if dat is null
  if (is.null(dat)) {
    warning("Table is empty")
    return()
  }
  dat %>%
    tibble::as_tibble() %>%
    mygo::attach_goterm_genecount() %>%
    dplyr::select(-pvalue, -qvalue) %>%
    # Move geneID to the very last row
    dplyr::select(dplyr::everything(), geneID) %>%
    rmyknife::dt_datatable() %>%
    return()
}

#' Create KEGG pathway as PNG
#' @param dat dataframe containing columns EntrezID and log2FoldChange
#' @param pathway kegg pathway name
#' @param species kegg species ID
#' @export
#' @import pathview magrittr dplyr png
#' @return png object
plot_kegg <- function(
  dat,
  pathway,
  species = "mmu"
) {
  # Create a named vector containing foldchanges and Entrez IDs
  fc_symbol_entrez <- dat %>%
    dplyr::filter(!is.na(EntrezID)) %>%
    .$log2FoldChange %>%
    as.double() %>%
    `names<-`(dat %>% .$EntrezID)
  # https://yulab-smu.top/clusterProfiler-book/chapter12.html
  kegg_pathway <- pathview::pathview(
    gene.data  = fc_symbol_entrez,
    pathway.id = pathway,
    gene.gene.idtype = "entrez",
    species = species,
    low = list("gene" = "blue"),
    high = list("gene" = "red")
  )
  # Filename
  kegg_png <- png::readPNG(glue::glue("{getwd()}/{pathway}.pathview.png"))
  # Remove intermediate files
  file.remove(glue::glue("{getwd()}/{pathway}.pathview.png"))
  file.remove(glue::glue("{getwd()}/{pathway}.png"))
  file.remove(glue::glue("{pathway}.pathview.png"))
  file.remove(glue::glue("{pathway}.png"))
  file.remove(glue::glue("{pathway}.xml"))
  # For some reason elusive to me we cannot print it in here or export it
  # As an raster image without weird stuff happening. That's why we made
  # the printing png function separate.
  kegg_png %>%
    return()
}

#' Plot png object (e.g. from plot_kegg_pathway)
#' @param png object
#' @import png magrittr
#' @export
plot_png <- function(png) {
  png %>%
    as.raster() %>%
    plot()
}
