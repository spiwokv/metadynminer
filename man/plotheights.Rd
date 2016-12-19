\name{plotheights}
\alias{plotheights}
\title{
  Plots the profile of hill heights metadynamics
}
\description{
  Plots hill heights (scaled by (Delta T + T)/Delta T) as a function of time.
}
\usage{
plotheights(hills)
}
\arguments{
  \item{hills}{
    hillsfile is a data frame containing time, CV(s), hill width(s), height and bias factor.
  }
}
\details{
   Other plot parameters such as \code{xlab}, \code{ylab}, \code{xlim}, \code{ylim}, \code{zlim}, \code{main}, \code{sub}, \code{pch}, \code{col}, \code{bg}, \code{cex}, \code{asp}, \code{lwd} or \code{axes} can be used.
}
\value{
}
\references{
}
\author{
  Vojtech Spiwok
}
\note{
}

\seealso{
}
\examples{
  data(acealanme1.3)
  plotheights(acealanme1.3)
}
