% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools.R
\name{perform_enrichGO}
\alias{perform_enrichGO}
\title{Perform a clusterPrfiler enrichGO analysis}
\usage{
perform_enrichGO(
  ontology,
  entrezgenes,
  background_genes,
  use_background,
  species
)
}
\arguments{
\item{ontology}{Can be either "BP", "CC", "MF"}

\item{entrezgenes}{List of entrezgenes to use for GO analysis}

\item{use_background}{use specified background genes}

\item{species}{Species to use for GO analysis, either "HUM" or "MUS"}

\item{entrez_background_genes}{List of background genes}
}
\description{
Perform a clusterPrfiler enrichGO analysis
}
