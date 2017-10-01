library(Rcpp)
library(RcppArmadillo)

# read HILLS from Plumed
read.hills<-function(file="HILLS", per=c(FALSE, FALSE)) {
  hillsf<-read.table(file, header=F, comment.char="#")
  if(ncol(hillsf)==5 || ncol(hillsf)==6) {
    cat("1D HILLS file read\n")
    hills<-list(hillsfile=hillsf, size=dim(hillsf), filename=file, per=per)
    class(hills) <- "hillsfile"
    return(hills)
  } else {
    if(ncol(hillsf)==7 || ncol(hillsf)==8) {
      cat("2D HILLS file read\n")
      hills<-list(hillsfile=hillsf, size=dim(hillsf), filename=file, per=per)
      class(hills) <- "hillsfile"
      return(hills)
    } else {
      stop("number of columns in HILLS file must be 5 or 6 (1D) or 7 or 8 (3D)")
    }
  }
}

# print a hillsfile
print.hillsfile<-function(hills) {
  if(hills$size[2]==5) {
    cat("1D hills file ")
    cat(hills$filename)
    cat(" with ")
    cat(hills$size[1])
    cat(" lines\n")
  }
  if(hills$size[2]==7) {
    cat("2D hills file ")
    cat(hills$filename)
    cat(" with ")
    cat(hills$size[1])
    cat(" lines\n")
  }
}

# print summary of a hillsfile
summary.hillsfile<-function(hills) {
  if(hills$size[2]==5) {
    cat("1D hills file ")
    cat(hills$filename)
    cat(" with ")
    cat(hills$size[1])
    cat(" lines\n")
    cat("The CV1 ranges from ")
    cat(min(hills$hillsfile[,2]))
    cat(" to ")
    cat(max(hills$hillsfile[,2]))
    cat("\n")
  }
  if(hills$size[2]==7) {
    cat("2D hills file ")
    cat(hills$filename)
    cat(" with ")
    cat(hills$size[1])
    cat(" lines\n")
    cat("The CV1 ranges from ")
    cat(min(hills$hillsfile[,2]))
    cat(" to ")
    cat(max(hills$hillsfile[,2]))
    cat("\nThe CV2 ranges from ")
    cat(min(hills$hillsfile[,3]))
    cat(" to ")
    cat(max(hills$hillsfile[,3]))
    cat("\n")
  }
}

# plot hillsfile
plot.hillsfile<-function(hills, perCV1r=c(-pi,pi), perCV2r=c(-pi,pi),
                         xlab=NULL, ylab=NULL,
                         xlim=NULL, ylim=NULL,
                         main=NULL, sub=NULL,
                         pch=1, col="black", bg="red", cex=1,
                         asp=NULL, lwd=1, axes=TRUE) {
  xlims<-NULL
  ylims<-NULL
  if(!is.null(xlim)) {xlims<-xlim}
  if((hills$per[1]==T)&is.null(xlim)) {xlims<-perCV1r}
  if(!is.null(ylim)) {ylims<-ylim}
  if((hills$per[2]==T)&is.null(ylim)) {ylims<-perCV2r}
  if(hills$size[2]==5) {
    if(is.null(xlab)) xlab="time"
    if(is.null(ylab)) ylab="CV"
    plot(hills$hillsfile[,1], hills$hillsfile[,2], type="l",
         xlab=xlab, ylab=ylab,
         main=main, sub=sub,
         xlim=xlims, ylim=ylims,
         col=col, cex=cex, lwd=lwd,
         asp=asp, axes=axes)
  }
  if(hills$size[2]==7) {
    if(is.null(xlab)) xlab="CV1"
    if(is.null(ylab)) ylab="CV2"
    plot(hills$hillsfile[,2], hills$hillsfile[,3], type="p",
         xlab=xlab, ylab=ylab,
         main=main, sub=sub,
         xlim=xlims, ylim=ylims,
         pch=pch, col=col, bg=bg, cex=cex, lwd=lwd,
         asp=asp, axes=axes)
  }
}

# plot heights
plotheights<-function(hills, xlab=NULL, ylab=NULL,
                      xlim=NULL, ylim=NULL, zlim=NULL,
                      main=NULL, sub=NULL,
                      pch=1, col="black", bg="red", cex=1,
                      asp=NULL, lwd=1, axes=TRUE) {
  if(class(hills)=="hillsfile") {
    if(is.null(xlab)) xlab="time"
    if(is.null(ylab)) ylab="hill height"
    if(hills$size[2]==5) {
      plot(hills$hillsfile[,1], hills$hillsfile[,4], type="l",
           xlab=xlab, ylab=ylab,
           main=main, sub=sub,
           col=col, cex=cex, lwd=lwd,
           asp=asp, axes=axes)
    }
    if(hills$size[2]==7) {
      plot(hills$hillsfile[,1], hills$hillsfile[,6], type="l",
           xlab=xlab, ylab=ylab,
           main=main, sub=sub,
           col=col, cex=cex, lwd=lwd,
           asp=asp, axes=axes)
    }
  } else {
    stop("function plotheights requires object hillsfile as an input")
  }
}

# red FES from MetadynView
read.fes<-function(filename, dimension=2, per=c(TRUE, TRUE)) {
  ifile<-read.table(filename)
  rows<-sqrt(nrow(ifile))
  fes<-matrix(ifile[,3], nrow=rows)
  fes<-max(fes)-fes
  x<-min(ifile[,1])+(max(ifile[,1])-min(ifile[,1]))*0:(rows-1)/(rows-1)
  y<-min(ifile[,2])+(max(ifile[,2])-min(ifile[,2]))*0:(rows-1)/(rows-1)
  cfes<-list(fes=fes, rows=rows, dimension=dimension, per=per, x=x, y=y)
  class(cfes) <- "fes"
  return(cfes)
}

# calculate 2d fes by bias sum algorithm
fes2d<-function(hills, perCV1r=c(-pi,pi), perCV2r=c(-pi,pi),
                tmin=0, tmax=NULL, xlim=NULL, ylim=NULL, npoints=256) {
  if(hills$size[2]==5) {
    stop("It looks like a 1D FES, use fes1d instead")
  }
  if(!is.null(tmax)) {
    if(hills$size[1]<tmax) {
      cat("You requested more hills by tmax than available, using all hills\n")
      tmax<-hills$size[1]
    }
  }
  if(is.null(tmax)) {
    tmax<-hills$size[1]
  }
  if(max(hills$hillsfile[,4])/min(hills$hillsfile[,4])>1.00000000001) {
    stop("Bias Sum algorithm works only with hills of the same sizes")
  }
  if(max(hills$hillsfile[,5])/min(hills$hillsfile[,5])>1.00000000001) {
    stop("Bias Sum algorithm works only with hills of the same sizes")
  }
  sourceCpp("../src/mm.cpp")
  minCV1 <- min(hills$hillsfile[,2])
  maxCV1 <- max(hills$hillsfile[,2])
  minCV2 <- min(hills$hillsfile[,3])
  maxCV2 <- max(hills$hillsfile[,3])
  xlims<-c(minCV1-0.05*(maxCV1-minCV1), maxCV1+0.05*(maxCV1-minCV1))
  ylims<-c(minCV2-0.05*(maxCV2-minCV2), maxCV2+0.05*(maxCV2-minCV2))
  if(!is.null(xlim)) {xlims<-xlim}
  if((hills$per[1]==T)&is.null(xlim)) {xlims<-perCV1r}
  if(!is.null(ylim)) {ylims<-ylim}
  if((hills$per[2]==T)&is.null(ylim)) {ylims<-perCV2r}
  x<-0:(npoints-1)*(xlims[2]-xlims[1])/(npoints-1)+xlims[1]
  y<-0:(npoints-1)*(ylims[2]-ylims[1])/(npoints-1)+ylims[1]
  if((hills$per[1]==F)&(hills$per[2]==F)) {
    fesm<-hills1(npoints*(hills$hillsfile[,2]-xlims[1])/(xlims[2]-xlims[1]),
                 npoints*(hills$hillsfile[,3]-ylims[1])/(ylims[2]-ylims[1]),
                 npoints*max(hills$hillsfile[,4])/(xlims[2]-xlims[1]),
                 npoints*max(hills$hillsfile[,5])/(ylims[2]-ylims[1]),
                 hills$hillsfile[,6],npoints,tmin,tmax)
  }
  if((hills$per[1]==T)&(hills$per[2]==F)) {
    fesm<-hills1p1(npoints*(hills$hillsfile[,2]-xlims[1])/(xlims[2]-xlims[1]),
                   npoints*(hills$hillsfile[,3]-ylims[1])/(ylims[2]-ylims[1]),
                   npoints*max(hills$hillsfile[,4])/(xlims[2]-xlims[1]),
                   npoints*max(hills$hillsfile[,5])/(ylims[2]-ylims[1]),
                   hills$hillsfile[,6],npoints,tmin,tmax)
  }
  if((hills$per[1]==F)&(hills$per[2]==T)) {
    fesm<-hills1p2(npoints*(hills$hillsfile[,2]-xlims[1])/(xlims[2]-xlims[1]),
                   npoints*(hills$hillsfile[,3]-ylims[1])/(ylims[2]-ylims[1]),
                   npoints*max(hills$hillsfile[,4])/(xlims[2]-xlims[1]),
                   npoints*max(hills$hillsfile[,5])/(ylims[2]-ylims[1]),
                   hills$hillsfile[,6],npoints,tmin,tmax)
  }
  if((hills$per[1]==T)&(hills$per[2]==T)) {
    fesm<-hills1p12(npoints*(hills$hillsfile[,2]-xlims[1])/(xlims[2]-xlims[1]),
                    npoints*(hills$hillsfile[,3]-ylims[1])/(ylims[2]-ylims[1]),
                    npoints*max(hills$hillsfile[,4])/(xlims[2]-xlims[1]),
                    npoints*max(hills$hillsfile[,5])/(ylims[2]-ylims[1]),
                    hills$hillsfile[,6],npoints,tmin,tmax)
  }
  cfes<-list(fes=fesm, rows=npoints, dimension=2, per=hills$per, x=x, y=y)
  class(cfes) <- "fes"
  return(cfes)
  return(cfes)
}

# calculate 2d fes conventionally (slow)
fes2d2<-function(hills, perCV1r=c(-pi,pi), perCV2r=c(-pi,pi),
                 tmin=0, tmax=NULL, xlim=NULL, ylim=NULL, npoints=256) {
  if(hills$size[2]==5) {
    stop("It looks like a 1D FES, use fes1d2 instead")
  }
  if(!is.null(tmax)) {
    if(hills$size[1]<tmax) {
      cat("You requested more hills by tmax than available, using all hills\n")
      tmax<-hills$size[1]
    }
  }
  if(is.null(tmax)) {
    tmax<-hills$size[1]
  }
  sourceCpp("../src/mm.cpp")
  minCV1 <- min(hills$hillsfile[,2])
  maxCV1 <- max(hills$hillsfile[,2])
  minCV2 <- min(hills$hillsfile[,3])
  maxCV2 <- max(hills$hillsfile[,3])
  xlims<-c(minCV1-0.05*(maxCV1-minCV1), maxCV1+0.05*(maxCV1-minCV1))
  ylims<-c(minCV2-0.05*(maxCV2-minCV2), maxCV2+0.05*(maxCV2-minCV2))
  if(!is.null(xlim)) {xlims<-xlim}
  if((hills$per[1]==T)&is.null(xlim)) {xlims<-perCV1r}
  if(!is.null(ylim)) {ylims<-ylim}
  if((hills$per[2]==T)&is.null(ylim)) {ylims<-perCV2r}
  x<-0:(npoints-1)*(xlims[2]-xlims[1])/(npoints-1)+xlims[1]
  y<-0:(npoints-1)*(ylims[2]-ylims[1])/(npoints-1)+ylims[1]
  if((hills$per[1]==F)&(hills$per[2]==F)) {
    fesm<-hills2(npoints*(hills$hillsfile[,2]-xlims[1])/(xlims[2]-xlims[1]),
                 npoints*(hills$hillsfile[,3]-ylims[1])/(ylims[2]-ylims[1]),
                 npoints*hills$hillsfile[,4]/(xlims[2]-xlims[1]),
                 npoints*hills$hillsfile[,5]/(ylims[2]-ylims[1]),
                 hills$hillsfile[,6],npoints,tmin,tmax)
  }
  if((hills$per[1]==T)&(hills$per[2]==F)) {
    fesm<-hills2p1(npoints*(hills$hillsfile[,2]-xlims[1])/(xlims[2]-xlims[1]),
                   npoints*(hills$hillsfile[,3]-ylims[1])/(ylims[2]-ylims[1]),
                   npoints*hills$hillsfile[,4]/(xlims[2]-xlims[1]),
                   npoints*hills$hillsfile[,5]/(ylims[2]-ylims[1]),
                   hills$hillsfile[,6],npoints,tmin,tmax)
  }
  if((hills$per[1]==F)&(hills$per[2]==T)) {
    fesm<-hills2p2(npoints*(hills$hillsfile[,2]-xlims[1])/(xlims[2]-xlims[1]),
                   npoints*(hills$hillsfile[,3]-ylims[1])/(ylims[2]-ylims[1]),
                   npoints*hills$hillsfile[,4]/(xlims[2]-xlims[1]),
                   npoints*hills$hillsfile[,5]/(ylims[2]-ylims[1]),
                   hills$hillsfile[,6],npoints,tmin,tmax)
  }
  if((hills$per[1]==T)&(hills$per[2]==T)) {
    fesm<-hills2p12(npoints*(hills$hillsfile[,2]-xlims[1])/(xlims[2]-xlims[1]),
                    npoints*(hills$hillsfile[,3]-ylims[1])/(ylims[2]-ylims[1]),
                    npoints*hills$hillsfile[,4]/(xlims[2]-xlims[1]),
                    npoints*hills$hillsfile[,5]/(ylims[2]-ylims[1]),
                    hills$hillsfile[,6],npoints,tmin,tmax)
  }
  cfes<-list(fes=fesm, rows=npoints, dimension=2, per=hills$per, x=x, y=y)
  class(cfes) <- "fes"
  return(cfes)
  return(cfes)
}

# calculate 1d fes by bias sum algorithm
fes1d<-function(hills, perCV1r=c(-pi,pi),
                tmin=0, tmax=NULL, xlim=NULL, npoints=256) {
  if(hills$size[2]==7) {
    stop("It looks like a 2D FES, use fes2d instead")
  }
  if(!is.null(tmax)) {
    if(hills$size[1]<tmax) {
      cat("You requested more hills by tmax than available, using all hills\n")
      tmax<-hills$size[1]
    }
  }
  if(is.null(tmax)) {
    tmax<-hills$size[1]
  }
  if(max(hills$hillsfile[,3])/min(hills$hillsfile[,3])>1.00000000001) {
    stop("Bias Sum algorithm works only with hills of the same sizes")
  }
  sourceCpp("../src/mm.cpp")
  minCV1 <- min(hills$hillsfile[,2])
  maxCV1 <- max(hills$hillsfile[,2])
  xlims<-c(minCV1-0.05*(maxCV1-minCV1), maxCV1+0.05*(maxCV1-minCV1))
  if(!is.null(xlim)) {xlims<-xlim}
  if((hills$per[1]==T)&is.null(xlim)) {xlims<-perCV1r}
  x<-0:(npoints-1)*(xlims[2]-xlims[1])/(npoints-1)+xlims[1]
  if(hills$per[1]==F) {
    fesm<-hills1d1(npoints*(hills$hillsfile[,2]-xlims[1])/(xlims[2]-xlims[1]),
                   npoints*max(hills$hillsfile[,3])/(xlims[2]-xlims[1]),
                   hills$hillsfile[,4],npoints,tmin,tmax)
  }
  if(hills$per[1]==T) {
    fesm<-hills1d1p(npoints*(hills$hillsfile[,2]-xlims[1])/(xlims[2]-xlims[1]),
                    npoints*max(hills$hillsfile[,3])/(xlims[2]-xlims[1]),
                    hills$hillsfile[,4],npoints,tmin,tmax)
  }
  cfes<-list(fes=fesm, rows=npoints, dimension=1, per=hills$per, x=x)
  class(cfes) <- "fes"
  return(cfes)
  return(cfes)
}

# calculate 1d fes conventionally (slow)
fes1d2<-function(hills, perCV1r=c(-pi,pi),
                 tmin=0, tmax=NULL, xlim=NULL, npoints=256) {
  if(hills$size[2]==7) {
    stop("It looks like a 2D FES, use fes2d2 instead")
  }
  if(!is.null(tmax)) {
    if (hills$size[1]<tmax) {
      cat("You requested more hills by tmax than available, using all hills\n")
      tmax<-hills$size[1]
    }
  }
  if(is.null(tmax)) {
    tmax<-hills$size[1]
  }
  sourceCpp("../src/mm.cpp")
  minCV1 <- min(hills$hillsfile[,2])
  maxCV1 <- max(hills$hillsfile[,2])
  xlims<-c(minCV1-0.05*(maxCV1-minCV1), maxCV1+0.05*(maxCV1-minCV1))
  if(!is.null(xlim)) {xlims<-xlim}
  if((hills$per[1]==T)&is.null(xlim)) {xlims<-perCV1r}
  x<-0:(npoints-1)*(xlims[2]-xlims[1])/(npoints-1)+xlims[1]
  if(hills$per[1]==F) {
    fesm<-hills1d1(npoints*(hills$hillsfile[,2]-xlims[1])/(xlims[2]-xlims[1]),
                   npoints*max(hills$hillsfile[,3])/(xlims[2]-xlims[1]),
                   hills$hillsfile[,4],npoints,tmin,tmax)
  }
  if(hills$per[1]==T) {
    fesm<-hills1d1p(npoints*(hills$hillsfile[,2]-xlims[1])/(xlims[2]-xlims[1]),
                    npoints*max(hills$hillsfile[,3])/(xlims[2]-xlims[1]),
                    hills$hillsfile[,4],npoints,tmin,tmax)
  }
  cfes<-list(fes=fesm, rows=npoints, dimension=1, per=hills$per, x=x)
  class(cfes) <- "fes"
  return(cfes)
  return(cfes)
}

# sum fesses
`+.fes`<-function(fes1, fes2) {
  if((class(fes1)=="fes")&(class(fes2)=="fes")) {
    if(fes1$rows!=fes2$rows) {
      stop("free energy surfaces have different numbers of points, exiting")
    }
    if(fes1$dimension!=fes2$dimension) {
      stop("free energy surfaces have different dimension, exiting")
    }
    if(sum(fes1$x!=fes2$x)>0) {
      stop("free energy surfaces have different CV1 axes, exiting")
    }
    if(fes1$dimension==2) {
      if(sum(fes1$y!=fes2$y)>0) {
        stop("free energy surfaces have different CV2 axes, exiting")
      }
    }
    if(fes1$dimension==1) {
      cfes<-list(fes=fes1$fes+fes2$fes, rows=fes1$rows, dimension=fes1$dimension, per=fes1$per, x=fes1$x)
    }
    if(fes1$dimension==2) {
      cfes<-list(fes=fes1$fes+fes2$fes, rows=fes1$rows, dimension=fes1$dimension, per=fes1$per, x=fes1$x, y=fes1$y)
    }
  } else if(class(fes1)=="fes") {
    if(fes1$dimension==1) {
      cfes<-list(fes=fes1$fes+fes2, rows=fes1$rows, dimension=fes1$dimension, per=fes1$per, x=fes1$x)
    }
    if(fes1$dimension==2) {
      cfes<-list(fes=fes1$fes+fes2, rows=fes1$rows, dimension=fes1$dimension, per=fes1$per, x=fes1$x, y=fes1$y)
    }
  } else if(class(fes2)=="fes") {
    if(fes2$dimension==1) {
      cfes<-list(fes=fes1+fes2$fes, rows=fes2$rows, dimension=fes2$dimension, per=fes2$per, x=fes2$x)
    }
    if(fes2$dimension==2) {
      cfes<-list(fes=fes1+fes2$fes, rows=fes2$rows, dimension=fes2$dimension, per=fes2$per, x=fes2$x, y=fes2$y)
    }
  }
  class(cfes) <- "fes"
  return(cfes)
  return(cfes)
}

# substract fesses
`-.fes`<-function(fes1, fes2) {
  if((class(fes1)=="fes")&(class(fes2)=="fes")) {
    if(fes1$rows!=fes2$rows) {
      stop("free energy surfaces have different numbers of points, exiting")
    }
    if(fes1$dimension!=fes2$dimension) {
      stop("free energy surfaces have different dimension, exiting")
    }
    if(sum(fes1$x!=fes2$x)>0) {
      stop("free energy surfaces have different CV1 axes, exiting")
    }
    if(fes1$dimension==2) {
      if(sum(fes1$y!=fes2$y)>0) {
        stop("free energy surfaces have different CV2 axes, exiting")
      }
    }
    if(fes1$dimension==1) {
      cfes<-list(fes=fes1$fes-fes2$fes, rows=fes1$rows, dimension=fes1$dimension, per=fes1$per, x=fes1$x)
    }
    if(fes1$dimension==2) {
      cfes<-list(fes=fes1$fes-fes2$fes, rows=fes1$rows, dimension=fes1$dimension, per=fes1$per, x=fes1$x, y=fes1$y)
    }
  } else if(class(fes1)=="fes") {
    if(fes1$dimension==1) {
      cfes<-list(fes=fes1$fes-fes2, rows=fes1$rows, dimension=fes1$dimension, per=fes1$per, x=fes1$x)
    }
    if(fes1$dimension==2) {
      cfes<-list(fes=fes1$fes-fes2, rows=fes1$rows, dimension=fes1$dimension, per=fes1$per, x=fes1$x, y=fes1$y)
    }
  } else if(class(fes2)=="fes") {
    if(fes2$dimension==1) {
      cfes<-list(fes=fes1-fes2$fes, rows=fes2$rows, dimension=fes2$dimension, per=fes2$per, x=fes2$x)
    }
    if(fes2$dimension==2) {
      cfes<-list(fes=fes1-fes2$fes, rows=fes2$rows, dimension=fes2$dimension, per=fes2$per, x=fes2$x, y=fes2$y)
    }
  }
  class(cfes) <- "fes"
  return(cfes)
  return(cfes)
}

# multiply a fes
`*.fes`<-function(fes1, fes2) {
  if((class(fes1)=="fes")&(class(fes2)=="fes")) {
    stop("you cannot multiply fes by fes")
  } else if(class(fes1)=="fes") {
    if(fes1$dimension==1) {
      cfes<-list(fes=fes1$fes*fes2, rows=fes1$rows, dimension=fes1$dimension, per=fes1$per, x=fes1$x)
    }
    if(fes1$dimension==2) {
      cfes<-list(fes=fes1$fes*fes2, rows=fes1$rows, dimension=fes1$dimension, per=fes1$per, x=fes1$x, y=fes1$y)
    }
  } else if(class(fes2)=="fes") {
    if(fes2$dimension==1) {
      cfes<-list(fes=fes1*fes2$fes, rows=fes2$rows, dimension=fes2$dimension, per=fes2$per, x=fes2$x)
    }
    if(fes2$dimension==2) {
      cfes<-list(fes=fes1*fes2$fes, rows=fes2$rows, dimension=fes2$dimension, per=fes2$per, x=fes2$x, y=fes2$y)
    }
  }
  class(cfes) <- "fes"
  return(cfes)
  return(cfes)
}

# divide a fes
`/.fes`<-function(fes1, coef) {
  if((class(fes1)=="fes")&(class(coef)=="fes")) {
    stop("you cannot divide fes by fes")
  } else if(class(fes1)=="fes") {
    if(fes1$dimension==1) {
      cfes<-list(fes=fes1$fes/coef, rows=fes1$rows, dimension=fes1$dimension, per=fes1$per, x=fes1$x)
    }
    if(fes1$dimension==2) {
      cfes<-list(fes=fes1$fes/coef, rows=fes1$rows, dimension=fes1$dimension, per=fes1$per, x=fes1$x, y=fes1$y)
    }
  } else if(class(coef)=="fes") {
    stop("you cannot divide something by fes")
  }
  class(cfes) <- "fes"
  return(cfes)
  return(cfes)
}

# min of fes
min.fes<-function(inputfes, na.rm=NULL) {
  return(min(inputfes$fes, na.rm=na.rm))
}

# max of fes
max.fes<-function(inputfes, na.rm=NULL) {
  return(max(inputfes$fes, na.rm=na.rm))
}

# mean of fes
mean.fes<-function(inputfes, na.rm=NULL) {
  return(mean(inputfes$fes, na.rm=na.rm))
}

# print FES
print.fes<-function(inputfes) {
  if(inputfes$dimension==1) {
    cat("1D free energy surface with ")
    cat(inputfes$rows)
    cat(" points, maximum ")
    cat(max(inputfes$fes))
    cat(" and minimum ")
    cat(min(inputfes$fes))
    cat("\n")
  }
  if(inputfes$dimension==2) {
    cat("2D free energy surface with ")
    cat(inputfes$rows)
    cat(" x ")
    cat(inputfes$rows)
    cat(" points, maximum ")
    cat(max(inputfes$fes))
    cat(" and minimum ")
    cat(min(inputfes$fes))
    cat("\n")
  }
}

# print summary of a FES
summary.fes<-function(inputfes) {
  if(inputfes$dimension==1) {
    cat("1D free energy surface with ")
    cat(inputfes$rows)
    cat(" points, maximum ")
    cat(max(inputfes$fes))
    cat(" and minimum ")
    cat(min(inputfes$fes))
    cat("\n")
  }
  if(inputfes$dimension==2) {
    cat("2D free energy surface with ")
    cat(inputfes$rows)
    cat(" x ")
    cat(inputfes$rows)
    cat(" points, maximum ")
    cat(max(inputfes$fes))
    cat(" and minimum ")
    cat(min(inputfes$fes))
    cat("\n")
  }
}

# plot FES
plot.fes<-function(inputfes, plottype="both",
                  x=NULL, y=NULL,
                  xlim=NULL, ylim=NULL, zlim=NULL,
                  main=NULL, sub=NULL,
                  xlab=NULL, ylab=NULL,
                  nlevels=10, levels=NULL,
                  col=rainbow(135)[100:1],
                  labels=NULL, labcex=0.6, drawlabels=TRUE,
                  method="flattest",
                  contcol=par("fg"), lty=par("lty"), lwd=par("lwd"),
                  axes=T) {
  fes<-inputfes$fes
  rows<-inputfes$rows
  if(inputfes$dimension==1) {
    if(is.null(x)) x<-inputfes$x
    if(is.null(xlab)) xlab="CV"
    if(is.null(ylab)) ylab="free energy"
    if(is.null(xlim)) xlim<-c(min(x),max(x))
    if(is.null(ylim)) {
      ylim<-range(pretty(range(fes)))
    }
    plot(x, fes, type="l",
        col=col, xlim=xlim, ylim=ylim,
        xlab=xlab, ylab=ylab, axes=axes,
        main=main, sub=sub)
  } else {
    if(is.null(x)) x<-inputfes$x
    if(is.null(y)) y<-inputfes$y
    if(is.null(xlab)) xlab="CV1"
    if(is.null(ylab)) ylab="CV2"
    if(is.null(zlim)) {
      zlim<-range(pretty(range(fes)))
    }
    if(is.null(levels)) {
      levels<-pretty(zlim, nlevels)
    }
    if(is.null(xlim)) xlim<-c(min(x),max(x))
    if(is.null(ylim)) ylim<-c(min(y),max(y))
    if(plottype=="image" || plottype=="both") {
      image(x, y, fes, zlim=zlim,
        col=col, xlim=xlim, ylim=ylim,
        xlab=xlab, ylab=ylab, axes=axes,
        main=main, sub=sub)
    }
    if(plottype=="contour") {
      contour(x, y, fes, zlim=zlim,
              nlevels=nlevels, levels=levels,
              labels=labels, labcex=labcex, drawlabels=drawlabels,
              method=method, col=contcol, lty=lty, lwd=lwd,
              main=main, sub=sub)
    }
    if(plottype=="both") {
      contour(x, y, fes, zlim=zlim,
              nlevels=nlevels, levels=levels,
              labels=labels, labcex=labcex, drawlabels=drawlabels,
              method=method, col=contcol, lty=lty, lwd=lwd, add=T)
    }
  }
}

# find minima of a FES
fesminima<-function(inputfes) {
  fes<-inputfes$fes
  rows<-inputfes$rows
  per<-inputfes$per
  minx<-c()
  miny<-c()
  for(i in 0:7) {
    ni<-i*16+0:17
    if(per[1]) {
      ni[ni==0]<-128
      ni[ni==129]<-1
    } else {
      ni<-ni[ni!=0]
      ni<-ni[ni!=129]
    }
    for(j in 0:7) {
      nj<-(j*16+0:17)
      if(per[2]) {
        nj<-nj[nj!=0]
        nj<-nj[nj!=129]
      } else {
        nj[nj==0]<-128
        nj[nj==129]<-1
      }
      binmin<-which(fes[ni,nj]==min(fes[ni,nj]), arr.ind = TRUE)
      if(binmin[1]!=1 && binmin[1]!=18 && binmin[2]!=1 && binmin[2]!=18) {
        minx<-c(minx,i*16+binmin[1]-1)
        miny<-c(miny,j*16+binmin[2]-1)
      }
    }
  }
  myLETTERS <- c(LETTERS, paste("A", LETTERS, sep=""), paste("B", LETTERS, sep=""))[1:length(minx)]
  minima<-data.frame(myLETTERS, minx, miny, inputfes$x[minx], inputfes$y[miny], fes[cbind(minx,miny)])
  names(minima) <- c("letter", "CV1bin", "CV2bin", "CV1", "CV2", "free_energy")
  minima <- minima[order(minima[,6]),]
  rownames(minima) <- seq(length=nrow(minima))
  minima[,1]<-myLETTERS
  minima<-list(minima=minima, fes=fes, rows=rows, dimension=inputfes$dimension, per=per, x=inputfes$x, y=inputfes$y)
  class(minima) <- "minima"
  return(minima)
}

# print minima of a FES
print.minima<-function(minims) {
  cat("$minima\n\n")
  print(minims$minima)
}

# print a summary of minima of a FES
summary.minima<-function(minims, temp=300, eunit="kJ/mol") {
  toprint <- minims$minima
  if(eunit=="kJ/mol") {
    toprint<-cbind(toprint, exp(-1000*toprint[,6]/8.314/temp))
  }
  if(eunit=="J/mol") {
    toprint<-cbind(toprint, exp(-toprint[,6]/8.314/temp))
  }
  if(eunit=="kcal/mol") {
    toprint<-cbind(toprint, exp(-1000*toprint[,6]/8.314/temp/4.184))
  }
  if(eunit=="cal/mol") {
    toprint<-cbind(toprint, exp(-toprint[,6]/8.314/temp/4.184))
  }
  sumpop<-sum(toprint[,7])
  toprint<-cbind(toprint, 100*toprint[,7]/sumpop)
  names(toprint)[7]<-"relative_pop"
  names(toprint)[8]<-"pop"
  print(toprint)
}

# plot minima
plot.minima <- function(minims, plottype="both",
                  x=NULL, y=NULL,
                  xlim=NULL, ylim=NULL, zlim=NULL,
                  main=NULL, sub=NULL,
                  xlab=NULL, ylab=NULL,
                  nlevels=10, levels=NULL,
                  col=rainbow(135)[100:1],
                  labels=NULL, labcex=0.6, drawlabels=TRUE,
                  method="flattest", textcol="black",
                  pch=1, bg="red", cex=1,
                  contcol=par("fg"), lty=par("lty"), lwd=par("lwd"),
                  axes=T) {
  fes<-minims$fes
  rows<-minims$rows
  minpoints<-minims$minima[,4:5]
  minlabs<-minims$minima[,1]
  if(minims$dimension==1) {
    if(is.null(x)) x<-minims$x
    if(is.null(xlab)) xlab="CV"
    if(is.null(ylab)) ylab="free energy"
    if(is.null(xlim)) xlim<-c(min(x),max(x))
    if(is.null(ylim)) {
      ylim<-range(pretty(range(fes)))
    }
    plot(x, fes, type="l",
        col=col, xlim=xlim, ylim=ylim,
        xlab=xlab, ylab=ylab, axes=axes,
        main=main, sub=sub)
  } else {
    if(is.null(x)) x<-minims$x
    if(is.null(y)) y<-minims$y
    if(is.null(xlab)) xlab="CV1"
    if(is.null(ylab)) ylab="CV2"
    if(is.null(zlim)) {
      zlim<-range(pretty(range(fes)))
    }
    if(is.null(levels)) {
      levels<-pretty(zlim, nlevels)
    }
    if(is.null(xlim)) xlim<-c(min(x),max(x))
    if(is.null(ylim)) ylim<-c(min(y),max(y))
    if(plottype=="points") {
      text(minpoints, labels=minlabs, col=textcol, xlim=xlim, ylim=ylim,
        xlab=xlab, ylab=ylab, axes=axes,
        pch=pch, bg=bg, cex=cex,
        main=main, sub=sub)
    }
    if(plottype=="image" || plottype=="both") {
      image(x, y, fes, zlim=zlim,
        col=col, xlim=xlim, ylim=ylim,
        xlab=xlab, ylab=ylab, axes=axes,
        main=main, sub=sub)
      text(minpoints, labels=minlabs, col=textcol, xlim=xlim, ylim=ylim,
        xlab=xlab, ylab=ylab,
        pch=pch, bg=bg, cex=cex,
        main=main, sub=sub)
    }
    if(plottype=="contour") {
      contour(x, y, fes, zlim=zlim,
              nlevels=nlevels, levels=levels,
              labels=labels, labcex=labcex, drawlabels=drawlabels,
              method=method, col=contcol, lty=lty, lwd=lwd,
              main=main, sub=sub)
      text(minpoints, labels=minlabs, col=textcol, xlim=xlim, ylim=ylim,
        xlab=xlab, ylab=ylab,
        pch=pch, bg=bg, cex=cex,
        main=main, sub=sub)
    }
    if(plottype=="both") {
      contour(x, y, fes, zlim=zlim,
              nlevels=nlevels, levels=levels,
              labels=labels, labcex=labcex, drawlabels=drawlabels,
              method=method, col=contcol, lty=lty, lwd=lwd, add=T)
    }
  }
}

