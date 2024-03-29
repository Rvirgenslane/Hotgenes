% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AutoTation.R
\name{AutoTation}
\alias{AutoTation}
\title{Builds color annotations for pheatmap}
\usage{
AutoTation(
  Design = NULL,
  VarIds = NULL,
  Colors = NULL,
  binary_var = c("white", "black")
)
}
\arguments{
\item{Design}{Hotgenes object or a data.frame containing experimental design data}

\item{VarIds}{Vector of Strings to select design variables. Default (NULL), returns}

\item{Colors}{vector for desired colors}

\item{binary_var}{vector for desired binary colors.}
}
\value{
List of colors matching design variables
}
\description{
Builds color annotations for pheatmap
}
\examples{
Example_Hotgenes_dir<-system.file("extdata",
"Example_Hotgenes.Rdata",
package = "Hotgenes", mustWork = TRUE)
load(Example_Hotgenes_dir)
AutoTation(Example_Hotgenes, 
Colors=c("cyan", "magenta", "blue", "pink"),
binary_var=c("white", "black"))
}
