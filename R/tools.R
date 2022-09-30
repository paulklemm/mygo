#' Make GO-Term analysis and print out Report
#'
#' @export
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
#' @param ontology Can be either "BP", "CC", "MF"
#' @param entrezgenes List of entrezgenes to use for GO analysis
#' @param entrez_background_genes List of background genes
#' @param use_background use specified background genes
#' @param species Species to use for GO analysis, either "HUM" or "MUS"
perform_enrichGO <- function(
  ontology,
  entrezgenes,
  background_genes,
  use_background,
  species
) {
  if (species == "MUS") {
    org_db <- org.Mm.eg.db::org.Mm.eg.db
  } else if (species == "HUM") {
    org_db <- org.Hs.eg.db::org.Hs.eg.db
  } else {
    paste0("Unknown species: ", species) %>%
      stop()
  }
  if (use_background) {
    clusterProfiler::enrichGO(
      gene = entrezgenes,
      OrgDb = org_db,
      universe = background_genes,
      ont = ontology,
      readable = TRUE
    ) %>%
    return()
  } else {
    clusterProfiler::enrichGO(
      gene = entrezgenes,
      OrgDb = org_db,
      ont = ontology,
      readable = TRUE
    ) %>%
    return()
  }
}

#' Return a volcano plot with annotated top_n values
#'
#' @export
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
#' @export
#' @param dat Data frame containing columns `pValue` and `EntrezID`
#' @param significance_cutoff Cutoff to consider genes as significant
#' @param use_background Use gene list as background instead of using all genes as background
#' @param species Either "HUM" or "MUS"
get_go_all_ontologies <- function(
  dat,
  significance_cutoff = 0.05,
  use_background = TRUE,
  species = "MUS"
) {
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
  
  get_go_all_ontologies_helper(
    significant_entrezgenes,
    background_entrezgenes,
    use_background,
    species
  ) %>%
    return()
}

#' Get GO-terms for all ontologies. This method requires a `significant` flag in the data frame
#'
#' @export
#' @param dat Data frame containing columns `significant` and `EntrezID`
#' @param use_background Use gene list as background instead of using all genes as background
#' @param species Either "HUM" or "MUS"
get_go_all_ontologies_2 <- function(
  dat,
  use_background = TRUE,
  species = "MUS"
) {
  # Prepare data frame
  valid_dat <-
    dat %>%
    dplyr::filter(!is.na(EntrezID)) %>%
    dplyr::mutate(EntrezID = as.character(EntrezID))

  # Get all genes in the data set as background universe
  background_entrezgenes <-
    valid_dat %>%
    .$EntrezID
  # Get significant genes
  significant_entrezgenes <-
    valid_dat %>%
    dplyr::filter(significant == TRUE) %>%
    .$EntrezID
  
  get_go_all_ontologies_helper(
    significant_entrezgenes,
    background_entrezgenes,
    use_background,
    species
  ) %>%
    return()
}

#' Helper function getting all GO-terms based on vectors of Entrez-IDs
#'
#' @export
#' @param significant_entrezgenes Vector of significant genes as Entrez-IDs
#' @param background_entrezgenes Vector of background genes as Entrez-IDs
#' @param use_background Use background entrezgenes Vector
#' @param species Either "HUM" or "MUS"
get_go_all_ontologies_helper <- function(
  significant_entrezgenes,
  background_entrezgenes,
  use_background,
  species
) {
  # Diagnostics Species output
  paste0("Get GO-terms for all ontologies for species ", species) %>%
    message()
  # Perform GO-term analysis for each ontology
  go_terms <- list()
  go_terms$Biological_Process <-
    significant_entrezgenes %>%
    perform_enrichGO(
      ontology = "BP",
      background_genes = background_entrezgenes,
      use_background = use_background,
      species = species
    )
  go_terms$Molecular_Function <-
    significant_entrezgenes %>%
    perform_enrichGO(
      ontology = "MF",
      background_genes = background_entrezgenes,
      use_background = use_background,
      species = species
    )
  go_terms$Cellular_Components <-
    significant_entrezgenes %>%
    perform_enrichGO(
      ontology = "CC",
      background_genes = background_entrezgenes,
      use_background = use_background,
      species = species
    )
  return(go_terms)
}

#' Print a couple of standard plots for provided GO ontologies
#'
#' @export
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
#' @export
#' @param dat Table containing columns `fc` and `EntrezID`
#' @param p_cutoff P-value cutoff for GSEA
#' @param set_readable Replace EntrezIDs by readable gene names. This causes problems with ridgeplot
#' @param simplify Simplify terms
#' @param simplify_cutoff cutoff value of clusterProfiler::simplify()
get_gse_all_ontologies <- function(
  dat,
  p_cutoff = 0.05,
  set_readable = TRUE,
  simplify = FALSE,
  simplify_cutoff = 0.7
) {
  # Prepare input data as required by GSE
  fc <- dat %>% get_named_fc_vector()
  gse_terms <- list()
  # Perform GSEA for all ontologies
  gse_terms$Biological_Process <- perform_gseGO(
    ontology = "BP",
    fc = fc,
    p_cutoff = p_cutoff,
    set_readable = set_readable
  )
  gse_terms$Molecular_Function <- perform_gseGO(
    ontology = "MF",
    fc = fc,
    p_cutoff = p_cutoff,
    set_readable = set_readable
  )
  gse_terms$Cellular_Components <- perform_gseGO(
    ontology = "CC",
    fc = fc,
    p_cutoff = p_cutoff,
    set_readable = set_readable
  )
  # Simplify if required
  if (simplify) {
    gse_terms <-
      gse_terms %>%
      purrr::map(function(dat) {
        clusterProfiler::simplify(dat, cutoff = simplify_cutoff)
      })
  }
  return(gse_terms)
}

#' Get KEGG GSE terms
#'
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
#' @param ontology Can be either "BP", "CC", "MF"
#' @param fc Named vector of foldchanges (name denotes Entrez ID)
#' @param p_cutoff P-value cutoff for GSEA
#' @param set_readable Replace EntrezIDs by readable gene names. This causes problems with ridgeplot
perform_gseGO <- function(ontology, fc, p_cutoff = 0.05, set_readable = TRUE) {
  gse <- clusterProfiler::gseGO(
    geneList = fc,
    OrgDb = org.Mm.eg.db,
    ont = ontology,
    pvalueCutoff = p_cutoff,
    verbose = FALSE
  )
  if (set_readable) {
    gse <-
      gse %>%
      DOSE::setReadable(OrgDb = org.Mm.eg.db)
  }
  return(gse)
}

#' Perform a clusterProfiler gseKEGG analysis
#'
#' @export
#' @param fc Named vector of foldchanges (name denotes Entrez ID)
perform_gseKEGG <- function(fc) {
  fc %>% clusterProfiler::gseKEGG(organism = "mmu") %>%
    return()
}

#' Plot total gene count against percentage of significant genes
#' @export
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
    dat <-
      dat %>%
      dplyr::mutate(Description = forcats::fct_reorder(Description, Percent_Significant))
    plot_title <- paste0("Top enriched terms ordered by highest overlap")
  } else if (order_by == "significance") {
    dat <-
      dat %>%
      dplyr::mutate(Description = forcats::fct_reorder(Description, dplyr::desc(p.adjust)))
    plot_title <- paste0("Top enriched terms ordered by significance")
  } else {
    warning(paste0("I don't know what to do with order_by parameter value ", order_by))
    return()
  }
  dat %>%
    dplyr::rename(`Significant Gene Count` = Count) %>%
    head(n) %>%
    ggplot2::ggplot(
      mapping = ggplot2::aes(
        Percent_Significant,
        Description
      )
    ) +
      ggplot2::geom_segment(ggplot2::aes(xend = 0, yend = Description)) +
      ggplot2::geom_point(ggplot2::aes(color = p.adjust, size = `Significant Gene Count`)) +
      ggplot2::scale_color_continuous(high = "#132B43", low = "#56B1F7") +
      ggplot2::scale_size_continuous(range = c(1, 5)) +
      ggplot2::xlab("Percentage of Significant Genes") +
      ggplot2::ylab(NULL) +
      ggplot2::ggtitle(plot_title)
}

#' Helper function, takes clusterProfiler GO-term result and prints it as DT widget
#' @export
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
#' @return png object
#' @examples
#'   mygo::plot_kegg(kegg_dat, pathway = "mmu04110") %>% mygo::plot_png()
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
#' @export
plot_png <- function(png) {
  png %>%
    as.raster() %>%
    plot()
}

#' Combine goterms into a single data frame
#' 
#' @export
#' @param dat clusterProfiler result
#' @return data frame with goterms and GO-term source column
bind_goterm_table <- function(dat) {
  dat %>%
    purrr::imap(function(domain, domain_name) {
      domain %>%
        tibble::as_tibble() %>%
        mygo::attach_goterm_genecount() %>%
        dplyr::select(dplyr::everything(), geneID) %>%
        dplyr::select(-c(pvalue, qvalue)) %>%
        dplyr::mutate(source = stringr::str_replace(domain_name, "_", " "))
    }) %>%
    dplyr::bind_rows()
}

#' Combine emap plots into one ggplot2 plot
#' 
#' @param dat clusterProfiler result
#' @param n See enrichplot::emapplot
#' @return ggplot2 plot
#' @export
emap_plot_facet_category <- function(dat, n = 50) {
  dat$Biological_Process %>%
    mygo::emap_plot("Enrich Map Biological Process", n) +
    dat$Cellular_Components %>%
    mygo::emap_plot("Enrich Map Cellular Component", n) +
    dat$Molecular_Function %>%
    mygo::emap_plot("Enrich Map Molecular Function", n)
}

#' Wrapper for mygo::overlap_percentage_plot to facet all GO categories
#' 
#' @param dat clusterProfiler result
#' @param order_by See mygo::overlap_percentage_plot
#' @param n See mygo::overlap_percentage_plot
#' @return ggplot2 object
#' @export
overlap_percentage_plot_facet_category <- function(
  dat,
  order_by = "significance",
  n = 20
) {
  dat %>%
    # Attach goterm gene count to all terms
    purrr::imap(~
      .x %>%
        tibble::as_tibble() %>%
        mygo::attach_goterm_genecount() %>%
        dplyr::mutate(source = .y)
    ) %>%
    # Bind into a single dataframe
    dplyr::bind_rows() %>%
    dplyr::group_by(source) %>%
    dplyr::slice_head(n = n) %>%
    mygo::overlap_percentage_plot(
      order_by = order_by,
      # n doesn't matter since we already sliced the data
      n = 1e10
    ) +
      ggplot2::facet_grid(source~., scales = "free")
}

#' Run GO-Term analysis and retrieve plots
#' This is for usage in custom analysis scripts
#' @param dat Table containing `ensembl_gene_id`, `padj` and optionally `log2FoldChange`
#' @param significance_cutoff Adjusted p-value cutoff
#' @param use_background Use whole table as gene background
#' @param simplify Summarise GO-terms
#' @param simplify_cutoff When summarising, use this cutoff
#' @param run_kegg Run KEGG pathway analysis. Adds columns `kegg` and `kegg_dt`
#' @param species Species string, either "MUS" or "HUM"
#' @return list containing plots and table outputs
#' @export
run_goterms <- function(
  dat,
  significance_cutoff = .05,
  use_background = FALSE,
  simplify = TRUE,
  simplify_cutoff = .7,
  run_kegg = FALSE,
  species = "MUS"
) {
  # Get memoised versions to speed up repeated analyses
  get_go_all_ontologies <- rmyknife::get_memoised(mygo::get_go_all_ontologies)
  simplify_ontologies <- rmyknife::get_memoised(mygo::simplify_ontologies)
  enrichKEGG <- rmyknife::get_memoised(clusterProfiler::enrichKEGG)

  # Check species
  if (species == "MUS") {
    species_kegg <- "mmu"
  } else if (species_kegg == "HUM") {
    species_kegg <- "hsa"
  } else {
    stop("Species must be either 'MUS' or 'HUM'")
  }

  result <- list()

  dat_entrez <-
    dat %>%
    dplyr::rename(q_value = padj) %>%
    rmyknife::ensembl_to_entrez(
      ensembl_id_name = "ensembl_gene_id",
      keep_only_rows_with_entrez = TRUE,
      drop_duplicates = TRUE
    ) %>%
    rmyknife::attach_gene_symbol_from_entrez()

  result$goterms <-
    dat_entrez %>%
    get_go_all_ontologies(
      use_background = use_background,
      significance_cutoff = significance_cutoff,
      species = species
    )
  
  # Optionally simplify the ontologies
  if (simplify) {
    result$goterms <-
      result$goterms %>%
      simplify_ontologies(cutoff = simplify_cutoff)
  }

  # Are the goterms valid?
  valid_goterms <- (result$goterms %>% length() > 0)
  if (!valid_goterms) {
    warning("We found no GO-terms.")
  }

  result$plot_emap_bp <-
    result$goterms$Biological_Process %>%
    mygo::emap_plot("Enrich Map Biological Process")
  
  result$plot_emap_cc <-
    result$goterms$Cellular_Components %>%
    mygo::emap_plot("Enrich Map Cellular Component")
  
  result$plot_emap_mf <-
    result$goterms$Molecular_Function %>%
    mygo::emap_plot("Enrich Map Molecular Function")
  
  result$dt <-
    result$goterms %>%
    mygo::bind_goterm_table() %>%
    rmyknife::dt_datatable()
  
  if (valid_goterms) {
    result$plot_percentage <-
      result$goterms %>%
      mygo::overlap_percentage_plot_facet_category(
        order_by = "significance",
        n = 20
      )
  }

  # Run KEGG pathway analysis
  if (run_kegg) {
    # Get Kegg Pathways
    result$kegg <-
      enrichKEGG(
        gene = dat_entrez$EntrezID,
        organism = species_kegg,
        pvalueCutoff = significance_cutoff
      )

    result$kegg_dt <-
      result$kegg %>%
      tibble::as_tibble() %>%
      rmyknife::dt_datatable(caption = "KEGG Pathways")
  }
  
  # Operations related to log2FoldChange
  has_log2fc <- !is.null(dat$log2FoldChange)
  if (has_log2fc & valid_goterms) {
    result$volcano <-
      dat %>%
      rmyknife::plot_volcano(
        min_padj = significance_cutoff,
        min_log2fc = 0
      )
    
    # Get genes per GO-term ordered by log2FoldChange
    result$goterm_genes <-
      result$goterms %>%
      mygo::bind_goterm_table() %>%
      dplyr::select(source, ID, Description, geneID) %>%
      dplyr::rename(external_gene_name = geneID) %>%
      tidyr::separate_rows(external_gene_name, sep = "/") %>%
      dplyr::left_join(
        dat %>%
          dplyr::select(ensembl_gene_id, external_gene_name, padj, log2FoldChange),
        by = "external_gene_name"
      ) %>%
      dplyr::group_by(ID) %>%
      dplyr::arrange(dplyr::desc(log2FoldChange), .by_group = TRUE) %>%
      # Remove grouping
      dplyr::ungroup()
  }

  return(result)
}

#' Get Table containing Kegg pathway information `pathway`, `ensembl_gene_id`, `kegg_name`
#' @param species Species string, either "MUS" or "HUM"
#' @return Table containing Kegg pathway information `pathway`, `ensembl_gene_id`, `kegg_name`
#' @export
#' @examples
#' get_kegg_table("MUS")
get_kegg_table <- function(species = "MUS") {
  # Get memoised versions
  keggGet <- rmyknife::get_memoised(KEGGREST::keggGet)
  keggLink <- rmyknife::get_memoised(KEGGREST::keggLink)
  # Get KEGG pathways for each species
  if (species == "MUS") {
    kegg_dat <-
      keggLink("pathway", "mmu") %>% 
      tibble::tibble(pathway = ., eg = sub("mmu:", "", names(.))) %>%
      dplyr::mutate(
        # symbol = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, eg, "SYMBOL", "ENTREZID"),
        ensembl_gene_id = AnnotationDbi::mapIds(org.Mm.eg.db::org.Mm.eg.db, eg, "ENSEMBL", "ENTREZID")
      )
  } else if (species == "HUM") {
    kegg_dat <-
      keggLink("pathway", "hsa") %>% 
      tibble::tibble(pathway = ., eg = sub("hsa:", "", names(.))) %>%
      dplyr::mutate(
        # symbol = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, eg, "SYMBOL", "ENTREZID"),
        ensembl_gene_id = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, eg, "ENSEMBL", "ENTREZID")
      )
  } else {
    paste0("Species ", species, " not supported") %>%
      stop()
  }

  # Get pathway names. Make sure we do not pull pathway names multiple times since its slow
  pathways <-
    kegg_dat %>%
    dplyr::distinct(pathway) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      kegg_name = keggGet(pathway)[[1]]$NAME
    )

  kegg_dat %>%
    dplyr::left_join(pathways, by = "pathway") %>%
    dplyr::rename(kegg_pathway = pathway)
}
