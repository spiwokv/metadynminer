\name{plot.hillsfile}
\alias{plot.hillsfile}
\title{
  Plots the profile of CVs of metadynamics
}
\description{
  Plots the profile of CVs of metadynamics. For metadynamics with one CV it plots CV as a function of time.
  For metadynamics with two CVs it plots CV2 as a function of CV1.
}
\usage{
plot.hillsfile(hills, xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL, zlim = NULL, main = NULL, sub = NULL, pch = 1, col = "black", bg = "red", cex = 1, asp = NULL, lwd = 1, axes = TRUE)
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
    if (hills$size[2] == 5) {
        if (is.null(xlab)) 
            xlab = "time"
        if (is.null(ylab)) 
            ylab = "CV"
        plot(hills$hillsfile[, 1], hills$hillsfile[, 2], type = "l", 
            xlab = xlab, ylab = ylab, main = main, sub = sub, 
            col = col, cex = cex, lwd = lwd, asp = asp, axes = axes)
    }
    if (hills$size[2] == 7) {
        if (is.null(xlab)) 
            xlab = "CV1"
        if (is.null(ylab)) 
            ylab = "CV2"
        plot(hills$hillsfile[, 2], hills$hillsfile[, 3], type = "p", 
            xlab = xlab, ylab = ylab, main = main, sub = sub, 
            pch = pch, col = col, bg = bg, cex = cex, lwd = lwd, 
            asp = asp, axes = axes)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
