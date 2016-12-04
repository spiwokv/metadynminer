\name{read.fes}
\alias{read.fes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
read.fes(filename, dimension = 2, per = c(TRUE, TRUE))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{
%%     ~~Describe \code{filename} here~~
}
  \item{dimension}{
%%     ~~Describe \code{dimension} here~~
}
  \item{per}{
%%     ~~Describe \code{per} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (filename, dimension = 2, per = c(TRUE, TRUE)) 
{
    ifile <- read.table(filename)
    rows <- sqrt(nrow(ifile))
    fes <- matrix(ifile[, 3], nrow = rows)
    fes <- max(fes) - fes
    x <- min(ifile[, 1]) + (max(ifile[, 1]) - min(ifile[, 1])) * 
        0:(rows - 1)/(rows - 1)
    y <- min(ifile[, 2]) + (max(ifile[, 2]) - min(ifile[, 2])) * 
        0:(rows - 1)/(rows - 1)
    cfes <- list(fes = fes, rows = rows, dimension = dimension, 
        per = per, x = x, y = y)
    class(cfes) <- "fes"
    return(cfes)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
