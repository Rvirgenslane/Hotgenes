% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Pheatmap_Viz.R
\name{Pheatmap_Viz}
\alias{Pheatmap_Viz}
\title{Provides expression data for select genes and samples}
\usage{
Pheatmap_Viz(
  DE_Input = NULL,
  hotList = NULL,
  samples_ids = 1,
  readouts = 1,
  ncut = NULL,
  selected_contrast = 1,
  lfc_cut = 0
)
}
\arguments{
\item{DE_Input}{Object generated by DEseq2_export or Limma_export functions.}

\item{hotList}{vector of gene symbols to select}

\item{samples_ids}{vector of sample ids to select}

\item{readouts}{index of normalized data to use.
See names(Hotgenes_input$Normalized_Expression).}

\item{ncut}{numeric value for the maximum number of genes to show.}

\item{selected_contrast}{name or numeric value to retrieve differential 
expression data. See names(Hotgenes_input$Output_DE).}

\item{lfc_cut}{absolute logfoldchange cutoff.}
}
\description{
Provides expression data for select genes and samples
}
