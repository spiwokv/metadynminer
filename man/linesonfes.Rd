% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nebs.R
\name{linesonfes}
\alias{linesonfes}
\title{Plot lines for Nudged Elastic Band projected onto free energy surface}
\usage{
linesonfes(nebpath = nebpath, col = "red", lwd = 1)
}
\arguments{
\item{nebpath}{nebpath object}
}
\description{
`points.nebpath` plots lines for free energy profile calculated by Nudged Elastic Band
projected onto free energy surface.
}
\examples{
tfes<-fes(acealanme)
minima<-fesminima(tfes)
nebAD<-neb(minima, min1="A", min2="D")
plot(minima)
linesonfes(nebAD)
}