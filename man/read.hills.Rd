\name{read.hills}
\alias{read.hills}
\title{HILLS reader}
\description{Function to read HILLS file produced by Plumed. It read HILLS file with one or two CVs and
             automatically recognizes its the number of CVs. The output is a hillsfile object.
}
\usage{
read.hills(file = "HILLS", per = c(FALSE, FALSE))
}
\arguments{
  \item{file}{
     the name of the file which the data are to be read from. 
     If it does not contain an absolute path, the file name is
     relative to the current working directory, \code{getwd()}.
     Compressed files are not supported.
    
     \code{file} can also be a complete URL. Default name is \code{HILLS}.
}
  \item{per}{
     describes periodicity of CV(s). It is a vector with one or two logical
     values. If for example both CVs are periodic from -pi to +pi it should
     be set to \code{c(TRUE, TRUE)}. It is \code{FALSE} by default.
}
}
\details{
  The file must be valid HILLS file from Plumed 2.
}
\value{
  An object \code{hillsfile} containing hill positions, heights and widths.
  
  An object of class \code{hillsfile} is a \code{list} containing
  following components:
  \item{hillsfile}{
    hillsfile is a data frame containing time, CV(s), hill width(s), height and bias factor.
  }
  \item{size}{
    size id a vector containing the number of lines and columns of hillsfile.
  }
  \item{filename}{
    filename is the name of the HILLS file.
  }
}
\references{
  http://www.plumed.org
}
\author{
  Vojtech Spiwok
}
\note{
}

\examples{
  myhills <- read.hills("HILLS", per=c(TRUE,TRUE))
  summary(myhills)
}
