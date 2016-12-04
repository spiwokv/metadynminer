\name{summary.minima}
\alias{summary.minima}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
summary.minima(minims, temp = 300, eunit = "kJ/mol")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{minims}{
%%     ~~Describe \code{minims} here~~
}
  \item{temp}{
%%     ~~Describe \code{temp} here~~
}
  \item{eunit}{
%%     ~~Describe \code{eunit} here~~
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
function (minims, temp = 300, eunit = "kJ/mol") 
{
    toprint <- minims$minima
    if (eunit == "kJ/mol") {
        toprint <- cbind(toprint, exp(-1000 * toprint[, 6]/8.314/temp))
    }
    if (eunit == "J/mol") {
        toprint <- cbind(toprint, exp(-toprint[, 6]/8.314/temp))
    }
    if (eunit == "kcal/mol") {
        toprint <- cbind(toprint, exp(-1000 * toprint[, 6]/8.314/temp/4.184))
    }
    if (eunit == "cal/mol") {
        toprint <- cbind(toprint, exp(-toprint[, 6]/8.314/temp/4.184))
    }
    sumpop <- sum(toprint[, 7])
    toprint <- cbind(toprint, 100 * toprint[, 7]/sumpop)
    names(toprint)[7] <- "relative_pop"
    names(toprint)[8] <- "pop"
    print(toprint)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
