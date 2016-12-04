\name{read.hills}
\alias{read.hills}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
read.hills(file = "HILLS", per = c(FALSE, FALSE))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
%%     ~~Describe \code{file} here~~
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
function (file = "HILLS", per = c(FALSE, FALSE)) 
{
    hillsf <- read.table(file, header = F, comment.char = "#")
    if (ncol(hillsf) == 5) {
        cat("1D HILLS file read\n")
        hills <- list(hillsfile = hillsf, size = dim(hillsf), 
            filename = file, per = per)
        class(hills) <- "hillsfile"
        return(hills)
    }
    else {
        if (ncol(hillsf) == 7) {
            cat("2D HILLS file read\n")
            hills <- list(hillsfile = hillsf, size = dim(hillsf), 
                filename = file, per = per)
            class(hills) <- "hillsfile"
            return(hills)
        }
        else {
            stop("number of columns in HILLS file must be 5 (1D) of 7 (3D)")
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
