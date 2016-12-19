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
%%     ~~Describe \code{hills} here~~
}
  \item{xlab}{
%%     ~~Describe \code{xlab} here~~
}
  \item{ylab}{
%%     ~~Describe \code{ylab} here~~
}
  \item{xlim}{
%%     ~~Describe \code{xlim} here~~
}
  \item{ylim}{
%%     ~~Describe \code{ylim} here~~
}
  \item{zlim}{
%%     ~~Describe \code{zlim} here~~
}
  \item{main}{
%%     ~~Describe \code{main} here~~
}
  \item{sub}{
%%     ~~Describe \code{sub} here~~
}
  \item{pch}{
%%     ~~Describe \code{pch} here~~
}
  \item{col}{
%%     ~~Describe \code{col} here~~
}
  \item{bg}{
%%     ~~Describe \code{bg} here~~
}
  \item{cex}{
%%     ~~Describe \code{cex} here~~
}
  \item{asp}{
%%     ~~Describe \code{asp} here~~
}
  \item{lwd}{
%%     ~~Describe \code{lwd} here~~
}
  \item{axes}{
%%     ~~Describe \code{axes} here~~
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
function (hills, xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL, 
    zlim = NULL, main = NULL, sub = NULL, pch = 1, col = "black", 
    bg = "red", cex = 1, asp = NULL, lwd = 1, axes = TRUE) 
{
    if (class(hills) == "hillsfile") {
        if (is.null(xlab)) 
            xlab = "time"
        if (is.null(ylab)) 
            ylab = "hill height"
        if (hills$size[2] == 5) {
            plot(hills$hillsfile[, 1], hills$hillsfile[, 4], 
                type = "l", xlab = xlab, ylab = ylab, main = main, 
                sub = sub, col = col, cex = cex, lwd = lwd, asp = asp, 
                axes = axes)
        }
        if (hills$size[2] == 7) {
            plot(hills$hillsfile[, 1], hills$hillsfile[, 6], 
                type = "l", xlab = xlab, ylab = ylab, main = main, 
                sub = sub, col = col, cex = cex, lwd = lwd, asp = asp, 
                axes = axes)
        }
    }
    else {
        stop("function plotheights requires object hillsfile as an input")
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
