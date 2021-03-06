% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/minimaandpaths.R
\name{fesminima.fes}
\alias{fesminima.fes}
\title{Find free energy minima in the fes object}
\usage{
\method{fesminima}{fes}(inputfes, nbins = 8)
}
\arguments{
\item{inputfes}{fes object.}

\item{nbins}{number of bins for each CV (default 8).}
}
\value{
minima object.
}
\description{
`fesminima.fes` finds free energy minima on 1D or 2D free energy surface.
The surface is divided by a 1D or 2D grid and minima are found for each
bin. Next the program determines whether the minimum of a bin is a local
minimum of the whole free energy surface. Free energy minima are labeled
constitutively by capital letters.
}
\examples{
tfes<-fes(acealanme, imax=5000)
minima<-fesminima(tfes)
minima
}
