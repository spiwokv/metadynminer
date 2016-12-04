\name{fesminima}
\alias{fesminima}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fesminima(inputfes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputfes}{
%%     ~~Describe \code{inputfes} here~~
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
function (inputfes) 
{
    fes <- inputfes$fes
    rows <- inputfes$rows
    per <- inputfes$per
    minx <- c()
    miny <- c()
    for (i in 0:7) {
        ni <- i * 16 + 0:17
        if (per[1]) {
            ni[ni == 0] <- 128
            ni[ni == 129] <- 1
        }
        else {
            ni <- ni[ni != 0]
            ni <- ni[ni != 129]
        }
        for (j in 0:7) {
            nj <- (j * 16 + 0:17)
            if (per[2]) {
                nj <- nj[nj != 0]
                nj <- nj[nj != 129]
            }
            else {
                nj[nj == 0] <- 128
                nj[nj == 129] <- 1
            }
            binmin <- which(fes[ni, nj] == min(fes[ni, nj]), 
                arr.ind = TRUE)
            if (binmin[1] != 1 && binmin[1] != 18 && binmin[2] != 
                1 && binmin[2] != 18) {
                minx <- c(minx, i * 16 + binmin[1] - 1)
                miny <- c(miny, j * 16 + binmin[2] - 1)
            }
        }
    }
    myLETTERS <- c(LETTERS, paste("A", LETTERS, sep = ""), paste("B", 
        LETTERS, sep = ""))[1:length(minx)]
    minima <- data.frame(myLETTERS, minx, miny, inputfes$x[minx], 
        inputfes$y[miny], fes[cbind(minx, miny)])
    names(minima) <- c("letter", "CV1bin", "CV2bin", "CV1", "CV2", 
        "free_energy")
    minima <- minima[order(minima[, 6]), ]
    rownames(minima) <- seq(length = nrow(minima))
    minima[, 1] <- myLETTERS
    minima <- list(minima = minima, fes = fes, rows = rows, dimension = inputfes$dimension, 
        per = per, x = inputfes$x, y = inputfes$y)
    class(minima) <- "minima"
    return(minima)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
