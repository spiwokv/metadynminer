#' Find transition path on free energy surface by Nudged Elastic Band method
#'
#' `neb` finds a transition path on free energy surface for a given pair of
#' minima. For a 1D surface it simply takes the free energy profile between the
#' two minima. For 2D surface it calculates the transition path by Nudged Elastic
#; Band method (https://doi.org/10.1063/1.1323224).
#'
#' @param minims minima object.
#' @param min1 starting minimum identifier (can be letter or
#'        index, default "A").
#' @param min2 final minimum identifier (can be letter or index,
#'        default "B").
#' @param nbins number of bins along Nudged Elastic Band (default 20).
#' @param nsteps number of Nudged Elastic Band iterations (default 100).
#' @param step Nudged Elastic Band iteration step (default 1).
#' @param k Nudged Elastic Band toughness (default 0.2).
#' @return NEB path
#'
#' @export
#' @examples
#' tfes<-fes(acealanme, imax=5000)
#' minima<-fesminima(tfes)
#' nebAD<-neb(minima, min1="A", min2="D", nsteps=20)
#' nebAD
neb<-function(minims=minims, min1="A", min2="B", nbins=20,
              nsteps=100, step=1.0, k=0.2) {
  fes<-minims$fes
  pcv1<-minims$pcv1
  pcv2<-minims$pcv2
  rows<-minims$rows
  myLETTERS <- c(LETTERS, paste("A", LETTERS, sep=""), paste("B", LETTERS, sep=""))[1:nrow(minims$minima)]
  if(min1 %in% myLETTERS) {
    min1 <- minims$minima[match(min1, myLETTERS),]
  } else if(min1 %in% 1:nrow(minims$minima)) {
    min1 <- minims$minima[min1,]
  }
  if(min2 %in% myLETTERS) {
    min2 <- minims$minima[match(min2, myLETTERS),]
  } else if(min2 %in% 1:nrow(minims$minima)) {
    min2 <- minims$minima[min2,]
  }
  if(minims$dimension==1) {
    stop("Error: Nudged Elastic Band is available only for 2D free energy surfaces, exiting")
  }
  if(minims$dimension==2) {
    align<-function(v1, v2) {
      dot <- v1[1]*v2[1]+v1[2]*v2[2]
      return(c(dot*v2[1], dot*v2[2]))
    }
    perp<-function(xm1, x, xp1) {
      d1 <- sqrt((x[,1]-xm1[,1])**2 + (x[,2]-xm1[,2])**2)
      d2 <- sqrt((xp1[,1]-x[,1])**2 + (xp1[,2]-x[,2])**2)
      d <- d1 + d2
      return(c((xp1[,1]-xm1[,1])/d, (xp1[,2]-xm1[,2])/d))
    }
    cosphi<-function(xm1, x, xp1) {
      d1 <- sqrt((x[,1]-xm1[,1])**2 + (x[,2]-xm1[,2])**2)
      d2 <- sqrt((xp1[,1]-x[,1])**2 + (xp1[,2]-x[,2])**2)
      return(((x[,1]-xm1[,1])*(xp1[,1]-x[,1])+(x[,2]-xm1[,2])*(xp1[,2]-x[,2]))/d1/d2)
    }
    force1<-function(fes, x, y) {
      fx<-(fes[x-1,y]-fes[x+1,y])/2.0
      fy<-(fes[x,y-1]-fes[x,y+1])/2.0
      return(c(fx,fy))
    }
    force2<-function(xm1, x, xp1, k) {
      fx <- k*(xm1[,1]+xp1[,1]-2.0*x[,1])
      fy <- k*(xm1[,2]+xp1[,2]-2.0*x[,2])
      return(c(fx,fy))
    }
    x1<-c(min1[1,2], min1[1,3])
    x2<-c(min2[1,2], min2[1,3])
    pathx <- x1[1]+0:nbins*(x2[1]-x1[1])/nbins
    pathy <- x1[2]+0:nbins*(x2[2]-x1[2])/nbins
    fespot <- c()
    for(i in 1:(nbins+1)) {
      fespot<-c(fespot,fes[pathx[i],pathy[i]])
    }
    path<-data.frame(pathx,pathy,fespot)
    newpath<-path
    for(i in 1:nsteps) {
      for(j in 1:(nbins-1)) {
         tau <- perp(path[j,], path[j+1,], path[j+2,])
         x <- round(path[j+1,1])
         y <- round(path[j+1,2])
         x <- min(x,nrow(fes)-1)
         y <- min(y,ncol(fes)-1)
         f1 <- force1(fes, x, y)
         f1a <- align(f1, tau)
         f1p <- f1-f1a
         f2 <- force2(path[j,],path[j+1,],path[j+2,],k)
         f2a <- align(f2, tau)
         f2p <- f2-f2a
         fphi <- 0.5*(1.0+cos(cosphi(path[j,], path[j+1,], path[j+2,])*pi))
         newpath[j+1,1] <- newpath[j+1,1] + step*(f1p[1]+f2a[1]+fphi*f2p[1])
         newpath[j+1,2] <- newpath[j+1,2] + step*(f1p[2]+f2a[2]+fphi*f2p[2])
         newpath[j+1,3] <- fes[max(min(newpath[j+1,1],nrow(fes)),1),max(min(newpath[j+1,2],ncol(fes)),1)]
      }
      path <- newpath
    }
    path[,1]<-path[,1]*(pcv1[2]-pcv1[1])/rows+pcv1[1]
    path[,2]<-path[,2]*(pcv2[2]-pcv2[1])/rows+pcv2[1]
  }
  cnebpath<-list(path=path, min1=min1, min2=min2)
  class(cnebpath) <- "nebpath"
  return(cnebpath)
}

#' Print Nudged Elastic Band minima
#'
#' `print.nebpath` prints the list minima for Nudged Elastic Band
#'
#' @param x nebpath object
#' @param ... further arguments passed to or from other methods.
#'
#' @export
#' @examples
#' tfes<-fes(acealanme, imax=5000)
#' minima<-fesminima(tfes)
#' nebAD<-neb(minima, min1="A", min2="D", nsteps=20)
#' nebAD
print.nebpath <- function(x,...) {
  cat("path between minima:\n")
  print(x$min1)
  print(x$min2)
}

#' Print summary for Nudged Elastic Band
#'
#' `print.nebpath` prints the list minima for Nudged Elastic Band, activation energies and
#' half lives calculated by Eyring equation (https://doi.org/10.1063/1.1749604).
#'
#' @param object nebpath object.
#' @param temp temperature in Kelvins.
#' @param eunit energy units (kJ/mol or kcal/mol, kJ/mol is default).
#' @param ... further arguments passed to or from other methods.
#'
#' @export
#' @examples
#' tfes<-fes(acealanme, imax=5000)
#' minima<-fesminima(tfes)
#' nebAD<-neb(minima, min1="A", min2="D", nsteps=20)
#' summary(nebAD)
summary.nebpath <- function(object, temp=300, eunit="kJ/mol",...) {
  nebpath<-object
  cat("path between minima:\n")
  print(nebpath$min1)
  print(nebpath$min2)
  direction <- c("->", "<-")
  deltag <- c(0,0)
  halflife <- c(0,0)
  units <- c(NULL,NULL)
  if(ncol(nebpath$path)==3) {
    activ <- max(nebpath$path[,3])-nebpath$path[1,3]
    deltag[1]<-activ
    rate <- 1.38064852E-023*temp*exp(-1000.0*activ/8.314459848/temp)/6.62607004E-034;
    halft <- log(2)/rate
    halflife[1] <- halft
    units[1] <- "s"
    if(halft < 1.0) {
      halflife[1] <- halft*1000.0
      units[1] <- "ms"
    }
    if(halft < 0.001) {
      halflife[1] <- halft*1000000.0
      units[1] <- "micros"
    }
    if(halft < 0.000001) {
      halflife[1] <- halft*1000000000.0
      units[1] <- "ns"
    }
    if(halft < 0.000000001) {
      halflife[1] <- halft*1000000000000.0
      units[1] <- "ps"
    }
    if(halft < 0.000000000001) {
      halflife[1] <- halft*1000000000000000.0
      units[1] <- "fs"
    }
    if(halft > 60.0) {
      halflife[1] <- halft/60.0
      units[1] <- "min"
    }
    if(halft > 3600.0) {
      halflife[1] <- halft/3600.0
      units[1] <- "hours"
    }
    if(halft > 86400.0) {
      halflife[1] <- halft/86400.0
      units[1] <- "days"
    }
    if(halft > 31557600.0) {
      halflife[1] <- halft/31557600.0
      units[1] <- "years"
    }
    activ <- max(nebpath$path[,3])-nebpath$path[nrow(nebpath$path),3]
    deltag[2]<-activ
    rate <- 1.38064852E-023*temp*exp(-1000.0*activ/8.314459848/temp)/6.62607004E-034;
    halft <- log(2)/rate
    halflife[2] <- halft
    units[2] <- "s"
    if(halft < 1.0) {
      halflife[2] <- halft*1000.0
      units[2] <- "ms"
    }
    if(halft < 0.001) {
      halflife[2] <- halft*1000000.0
      units[2] <- "micros"
    }
    if(halft < 0.000001) {
      halflife[2] <- halft*1000000000.0
      units[2] <- "ns"
    }
    if(halft < 0.000000001) {
      halflife[2] <- halft*1000000000000.0
      units[2] <- "ps"
    }
    if(halft < 0.000000000001) {
      halflife[2] <- halft*1000000000000000.0
      units[2] <- "fs"
    }
    if(halft > 60.0) {
      halflife[2] <- halft/60.0
      units[2] <- "min"
    }
    if(halft > 3600.0) {
      halflife[2] <- halft/3600.0
      units[2] <- "hours"
    }
    if(halft > 86400.0) {
      halflife[2] <- halft/86400.0
      units[2] <- "days"
    }
    if(halft > 31557600.0) {
      halflife[2] <- halft/31557600.0
      units[2] <- "years"
    }
  }
  if(ncol(nebpath$path)==2) {
    activ <- max(nebpath$path[,2])-nebpath$path[1,2]
    deltag[1]<-activ
    rate <- 1.38064852E-023*temp*exp(-1000.0*activ/8.314459848/temp)/6.62607004E-034;
    halft <- log(2)/rate
    halflife[1] <- halft
    units[1] <- "s"
    if(halft < 1.0) {
      halflife[1] <- halft*1000.0
      units[1] <- "ms"
    }
    if(halft < 0.001) {
      halflife[1] <- halft*1000000.0
      units[1] <- "micros"
    }
    if(halft < 0.000001) {
      halflife[1] <- halft*1000000000.0
      units[1] <- "ns"
    }
    if(halft < 0.000000001) {
      halflife[1] <- halft*1000000000000.0
      units[1] <- "ps"
    }
    if(halft < 0.000000000001) {
      halflife[1] <- halft*1000000000000000.0
      units[1] <- "fs"
    }
    if(halft > 60.0) {
      halflife[1] <- halft/60.0
      units[1] <- "min"
    }
    if(halft > 3600.0) {
      halflife[1] <- halft/3600.0
      units[1] <- "hours"
    }
    if(halft > 86400.0) {
      halflife[1] <- halft/86400.0
      units[1] <- "days"
    }
    if(halft > 31557600.0) {
      halflife[1] <- halft/31557600.0
      units[1] <- "years"
    }
    activ <- max(nebpath$path[,2])-nebpath$path[nrow(nebpath$path),2]
    deltag[2]<-activ
    rate <- 1.38064852E-023*temp*exp(-1000.0*activ/8.314459848/temp)/6.62607004E-034;
    halft <- log(2)/rate
    halflife[2] <- halft
    units[2] <- "s"
    if(halft < 1.0) {
      halflife[2] <- halft*1000.0
      units[2] <- "ms"
    }
    if(halft < 0.001) {
      halflife[2] <- halft*1000000.0
      units[2] <- "micros"
    }
    if(halft < 0.000001) {
      halflife[2] <- halft*1000000000.0
      units[2] <- "ns"
    }
    if(halft < 0.000000001) {
      halflife[2] <- halft*1000000000000.0
      units[2] <- "ps"
    }
    if(halft < 0.000000000001) {
      halflife[2] <- halft*1000000000000000.0
      units[2] <- "fs"
    }
    if(halft > 60.0) {
      halflife[2] <- halft/60.0
      units[2] <- "min"
    }
    if(halft > 3600.0) {
      halflife[2] <- halft/3600.0
      units[2] <- "hours"
    }
    if(halft > 86400.0) {
      halflife[2] <- halft/86400.0
      units[2] <- "days"
    }
    if(halft > 31557600.0) {
      halflife[2] <- halft/31557600.0
      units[2] <- "years"
    }
  }
  cat("with kinetics\n")
  rates<-data.frame(direction, deltag, halflife, units)
  return(rates)
}

#' Plot Nudged Elastic Band
#'
#' `plot.nebpath` plots free energy profile calculated by Nudged Elastic Band.
#'
#' @param x nebpath object.
#' @param main an overall title for the plot: see 'title'.
#' @param sub a sub title for the plot: see 'title'.
#' @param xlab a title for the x axis: see 'title'.
#' @param ylab a title for the y axis: see 'title'.
#' @param asp the y/x aspect ratio, see 'plot.window'.
#' @param col color code or name, see 'par'.
#' @param cex text expansion.
#' @param lwd line width for drawing symbols see 'par'.
#' @param xlim numeric vector of length 2, giving the x coordinates range.
#' @param ylim numeric vector of length 2, giving the y coordinates range.
#' @param axes a logical value indicating whether both axes should be drawn
#'        on the plot.
#' @param ... further arguments passed to or from other methods.
#'
#' @export
#' @examples
#' tfes<-fes(acealanme, imax=5000)
#' minima<-fesminima(tfes)
#' nebAD<-neb(minima, min1="A", min2="D", nsteps=20)
#' plot(nebAD)
plot.nebpath <- function(x,
                         xlim=NULL, ylim=NULL,
                         main=NULL, sub=NULL,
                         xlab="bin", ylab="free energy",
                         col="red", lwd=1, asp=NULL, cex=1,
                         axes=T,...) {
  nebpath<-x
  if(ncol(nebpath$path)==3) {
    if(is.null(ylim)) {
      ylim<-c(min(nebpath$path[,3]), min(nebpath$path[,3])+1.1*(max(nebpath$path[,3])-min(nebpath$path[,3])))
    }
    plot(nebpath$path[,3], type="l",
         xlim=xlim, ylim=ylim,
         main=main, sub=sub,
         xlab=xlab, ylab=ylab,
         col=col, lwd=lwd,
         axes=T, asp=asp)
    text(c(1), c(nebpath$path[1,3])+0.05*(max(nebpath$path[,3])-min(nebpath$path[,3])), labels=c(nebpath$min1[1]), cex=cex)
    text(c(nrow(nebpath$path)), c(nebpath$path[nrow(nebpath$path),3])+0.05*(max(nebpath$path[,3])-min(nebpath$path[,3])),
         labels=c(nebpath$min2[1]), cex=cex)
    text(c(which.max(nebpath$path[,3])), c(nebpath$path[which.max(nebpath$path[,3]),3])+0.05*(max(nebpath$path[,3])-min(nebpath$path[,3])),
         labels=c("#"), cex=cex)
  }
  if(ncol(nebpath$path)==2) {
    if(is.null(ylim)) {
      ylim<-c(min(nebpath$path[,2]), min(nebpath$path[,2])+1.1*(max(nebpath$path[,2])-min(nebpath$path[,2])))
    }
    plot(nebpath$path[,2], type="l",
         xlim=xlim, ylim=ylim,
         main=main, sub=sub,
         xlab=xlab, ylab=ylab,
         col=col, lwd=lwd,
         axes=T, asp=asp)
    text(c(1), c(nebpath$path[1,2])+0.05*(max(nebpath$path[,2])-min(nebpath$path[,2])), labels=c(nebpath$min1[1]), cex=cex)
    text(c(nrow(nebpath$path)), c(nebpath$path[nrow(nebpath$path),2])+0.05*(max(nebpath$path[,2])-min(nebpath$path[,2])),
         labels=c(nebpath$min2[1]), cex=cex)
    text(c(which.max(nebpath$path[,2])), c(nebpath$path[which.max(nebpath$path[,2]),2])+0.05*(max(nebpath$path[,2])-min(nebpath$path[,2])),
         labels=c("#"), cex=cex)
  }
}

#' Plot points for Nudged Elastic Band
#'
#' `points.nebpath` plots points for free energy profile calculated by Nudged Elastic Band.
#'
#' @param x nebpath object.
#' @param pch plotting 'character', i.e., symbol to use. See 'points'.
#' @param col color code or name, see 'par'.
#' @param bg background (fill) color for the open plot symbols given by
#'        'pch = 21:25'.
#' @param cex character (or symbol) expansion: a numerical vector. This
#'        works as a multiple of 'par("cex")'.
#' @param lwd line width for drawing symbols see 'par'.
#' @param ... further arguments passed to or from other methods.
#'
#' @export
#' @examples
#' tfes<-fes(acealanme, imax=5000)
#' minima<-fesminima(tfes)
#' nebAD<-neb(minima, min1="A", min2="D", nsteps=20)
#' plot(nebAD)
#' points(nebAD)
points.nebpath <- function(x, pch=NULL, cex=1, bg=NULL,
                           col="red", lwd=1,...) {
  nebpath<-x
  if(ncol(nebpath$path)==3) {
    points(nebpath$path[,3],
           pch=NULL, cex=1, bg=NULL,
           col=col, lwd=lwd)
  }
  if(ncol(nebpath$path)==2) {
    points(nebpath$path[,2], type="l",
           pch=NULL, cex=1, bg=NULL,
           col=col, lwd=lwd)
  }
}

#' Plot lines for Nudged Elastic Band
#'
#' `lines.nebpath` plots lines for free energy profile calculated by Nudged Elastic Band.
#'
#' @param x nebpath object.
#' @param col color code or name, see 'par'.
#' @param lwd line width for drawing symbols see 'par'.
#' @param ... further arguments passed to or from other methods.
#'
#' @export
#' @examples
#' tfes<-fes(acealanme, imax=5000)
#' minima<-fesminima(tfes)
#' nebAD<-neb(minima, min1="A", min2="D", nsteps=20)
#' plot(nebAD)
#' lines(nebAD, lwd=4)
lines.nebpath <- function(x,
                          col="red", lwd=1,...) {
  nebpath<-x
  if(ncol(nebpath$path)==3) {
    lines(nebpath$path[,3],
          col=col, lwd=lwd)
  }
  if(ncol(nebpath$path)==2) {
    lines(nebpath$path[,2], type="l",
          col=col, lwd=lwd)
  }
}

#' Plot points for Nudged Elastic Band projected onto free energy surface
#'
#' `pointsonfes` plots points for free energy profile calculated by Nudged Elastic Band
#' projected onto free energy surface.
#'
#' @param x nebpath object.
#' @param pch plotting 'character', i.e., symbol to use. See 'points'.
#' @param col color code or name, see 'par'.
#' @param bg background (fill) color for the open plot symbols given by
#'        'pch = 21:25'.
#' @param cex character (or symbol) expansion: a numerical vector. This
#'        works as a multiple of 'par("cex")'.
#' @param lwd line width for drawing symbols see 'par'.
#'
#' @export
#' @examples
#' tfes<-fes(acealanme, imax=5000)
#' minima<-fesminima(tfes)
#' nebAD<-neb(minima, min1="A", min2="D", nsteps=20)
#' plot(minima)
#' pointsonfes(nebAD)
pointsonfes <- function(x,
                        pch=NULL, cex=1, bg=NULL,
                        col="red", lwd=1) {
  nebpath<-x
  if(ncol(nebpath$path)==3) {
    points(nebpath$path[,1], nebpath$path[,2],
           pch=NULL, cex=1, bg=NULL,
           col=col, lwd=lwd)
  }
  if(ncol(nebpath$path)==2) {
    stop("Error: Can be used only on 2D free energy surfaces")
  }
}

#' Plot lines for Nudged Elastic Band projected onto free energy surface
#'
#' `linesonfes` plots lines for free energy profile calculated by Nudged Elastic Band
#' projected onto free energy surface.
#'
#' @param x nebpath object.
#' @param col color code or name, see 'par'.
#' @param lwd line width for drawing symbols see 'par'.
#'
#' @export
#' @examples
#' tfes<-fes(acealanme, imax=5000)
#' minima<-fesminima(tfes)
#' nebAD<-neb(minima, min1="A", min2="D", nsteps=20)
#' plot(minima)
#' linesonfes(nebAD)
linesonfes <- function(x,
                       col="red", lwd=1) {
  nebpath<-x
  if(ncol(nebpath$path)==3) {
    lines(nebpath$path[,1], nebpath$path[,2],
          col=col, lwd=lwd)
  }
  if(ncol(nebpath$path)==2) {
    stop("Error: Can be used only on 2D free energy surfaces")
  }
}

