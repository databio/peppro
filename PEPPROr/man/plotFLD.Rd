% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PEPPROr.R
\name{plotFLD}
\alias{plotFLD}
\title{Plot fragment length distribution}
\usage{
plotFLD(
  fragL,
  fragL_count,
  fragL_txt = "fragLenDistribution.txt",
  max_fragment = 200
)
}
\arguments{
\item{fragL}{infile containing single column of fragment lengths}

\item{fragL_count}{counts of each fragment length identified}

\item{fragL_dis1}{pdf filename}

\item{ragL_dis2}{fragment length distribution stats file}
}
\description{
This function plots the fragment length distribution of a paired-end sample
and produces pdf/png files.
}
\examples{
data("frag_len")
data("frag_len_count")
plotFLD(fragL = "frag_len", fragL_count = "frag_len_count",
        fragL_txt = "fragLenDistribution_example.txt", max_fragment=200)
}
\keyword{distribution}
\keyword{fragment}
