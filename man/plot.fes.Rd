\name{plot.fes}
\alias{plot.fes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.fes(inputfes, plottype = "both", xlim = NULL, ylim = NULL, zlim = NULL, main = NULL, sub = NULL, xlab = NULL, ylab = NULL, nlevels = 10, levels = NULL, col = rainbow(135)[100:1], labels = NULL, labcex = 0.6, drawlabels = TRUE, method = "flattest", contcol = par("fg"), lty = par("lty"), lwd = par("lwd"), axes = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputfes}{
%%     ~~Describe \code{inputfes} here~~
}
  \item{plottype}{
%%     ~~Describe \code{plottype} here~~
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
  \item{xlab}{
%%     ~~Describe \code{xlab} here~~
}
  \item{ylab}{
%%     ~~Describe \code{ylab} here~~
}
  \item{nlevels}{
%%     ~~Describe \code{nlevels} here~~
}
  \item{levels}{
%%     ~~Describe \code{levels} here~~
}
  \item{col}{
%%     ~~Describe \code{col} here~~
}
  \item{labels}{
%%     ~~Describe \code{labels} here~~
}
  \item{labcex}{
%%     ~~Describe \code{labcex} here~~
}
  \item{drawlabels}{
%%     ~~Describe \code{drawlabels} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{contcol}{
%%     ~~Describe \code{contcol} here~~
}
  \item{lty}{
%%     ~~Describe \code{lty} here~~
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
function (inputfes, plottype = "both", xlim = NULL, ylim = NULL, 
    zlim = NULL, main = NULL, sub = NULL, xlab = NULL, ylab = NULL, 
    nlevels = 10, levels = NULL, col = rainbow(135)[100:1], labels = NULL, 
    labcex = 0.6, drawlabels = TRUE, method = "flattest", contcol = par("fg"), 
    lty = par("lty"), lwd = par("lwd"), axes = T) 
{
    fes <- inputfes$fes
    rows <- inputfes$rows
    if (inputfes$dimension == 1) {
        x <- inputfes$x
        if (is.null(xlab)) 
            xlab = "CV"
        if (is.null(ylab)) 
            ylab = "free energy"
        if (is.null(xlim)) 
            xlim <- c(min(x), max(x))
        if (is.null(ylim)) {
            ylim <- range(pretty(range(fes)))
        }
        plot(x, fes, type = "l", col = col, xlim = xlim, ylim = ylim, 
            xlab = xlab, ylab = ylab, axes = axes, main = main, 
            sub = sub)
    }
    else {
        x <- inputfes$x
        y <- inputfes$y
        if (is.null(xlab)) 
            xlab = "CV1"
        if (is.null(ylab)) 
            ylab = "CV2"
        if (is.null(zlim)) {
            zlim <- range(pretty(range(fes)))
        }
        if (is.null(levels)) {
            levels <- pretty(zlim, nlevels)
        }
        if (is.null(xlim)) 
            xlim <- c(min(x), max(x))
        if (is.null(ylim)) 
            ylim <- c(min(y), max(y))
        if (plottype == "image" || plottype == "both") {
            image(x, y, fes, zlim = zlim, col = col, xlim = xlim, 
                ylim = ylim, xlab = xlab, ylab = ylab, axes = axes, 
                main = main, sub = sub)
        }
        if (plottype == "contour") {
            contour(x, y, fes, zlim = zlim, nlevels = nlevels, 
                levels = levels, labels = labels, labcex = labcex, 
                drawlabels = drawlabels, method = method, col = contcol, 
                lty = lty, lwd = lwd, main = main, sub = sub)
        }
        if (plottype == "both") {
            contour(x, y, fes, zlim = zlim, nlevels = nlevels, 
                levels = levels, labels = labels, labcex = labcex, 
                drawlabels = drawlabels, method = method, col = contcol, 
                lty = lty, lwd = lwd, add = T)
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
