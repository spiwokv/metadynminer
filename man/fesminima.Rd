\name{fesminima}
\alias{fesminima}
\title{FES minima locator}
}
\description{Function to find minima on free energy surfaces calculated by metadynamics.
             The free energy surface is divided by 8x8 grid and minimum is located for each section.
             The output is a \code{minima} object.}
\usage{
fesminima(inputfes)
}
\arguments{
  \item{inputfes}{
      The object \code{inputfes} is a free energy surface calculated by function \code{mtdfes}, \code{mtdfes2}
      or read by \code{read.fes}.
  }
}
\details{
  The input object must be valid type \code{inputfes}.
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

\examples{
  myhills <- read.hills("http://www.metadynamics.cz/metadynminer/data/HILLS", per=c(TRUE,TRUE))
  myfes<-mtdfes(myhills)
  myminima<-fesminima(myfes)
  summary(myminima)
}
