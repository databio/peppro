% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PEPPROr.R
\name{mRNAcontamination}
\alias{mRNAcontamination}
\title{Plot the distribution of genic exonRPKM/intronRPKM ratios}
\usage{
mRNAcontamination(
  rpkm,
  name = "mRNA contamination ratios",
  raw = FALSE,
  type = c("histogram", "boxplot", "violin"),
  annotate = TRUE
)
}
\arguments{
\item{rpkm}{A three column TSV format file containing
"gene", "intron RPKM", "exon RPKM" columns.}

\item{name}{X-axis label, typically a sample name}

\item{raw}{Plot raw distribution}

\item{type}{Plot format}

\item{annotate}{Display mean and median values on plot}
}
\description{
This function plots the distribution of by gene exon RPKM divided by
intron RPKM ratios. Can produce raw or log10 distributions, but reports
both median values.
}
\examples{
data("rpkm_ratios")
mRNAcontamination(rpkm = "rpkm_ratios")
}
\keyword{contamination}
\keyword{mRNA}
