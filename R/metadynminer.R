library(Rcpp)
library(RcppArmadillo)

# read HILLS from Plumed
read.hills<-function(file="HILLS", per=c(FALSE, FALSE), pcv1=c(-pi,pi), pcv2=c(-pi,pi), ignoretime=FALSE) {
  hillsf<-read.table(file, header=F, comment.char="#")
  if(ncol(hillsf)==5 || ncol(hillsf)==6) {
    cat("1D HILLS file read\n")
    if(ignoretime) {
      if (hillsf[,1]!=sort(hillsf[,1])) {
        cat("Warning: time in the hills file is not continuous, probably you\n")
        cat("used RESTART function. The time will be updated automatically from zero\n")
        cat("according to the first step!\n")
        hillsf[,1]<-seq(from=hillsf[1,1], by=hillsf[1,1], length.out=nrow(hillsf))
      }
    }
    hills<-list(hillsfile=hillsf, size=dim(hillsf), filename=file, per=per)
    class(hills) <- "hillsfile"
    return(hills)
  } else {
    if(ncol(hillsf)==7 || ncol(hillsf)==8) {
      cat("2D HILLS file read\n")
      if(ignoretime) {
        if (hillsf[,1]!=sort(hillsf[,1])) {
          cat("Warning: time in the hills file is not continuous, probably you\n")
          cat("used RESTART function. The time will be updated automatically from zero\n")
          cat("according to the first step!\n")
          hillsf[,1]<-seq(from=hillsf[1,1], by=hillsf[1,1], length.out=nrow(hillsf))
        }
      }
      hills<-list(hillsfile=hillsf, size=dim(hillsf), filename=file, per=per, pcv1=pcv1, pcv2=pcv2)
      class(hills) <- "hillsfile"
      return(hills)
    } else {
      stop("number of columns in HILLS file must be 5 or 6 (1D) or 7 or 8 (2D)")
    }
  }
}

# print a hillsfile
print.hillsfile<-function(hills=hills) {
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
summary.hillsfile<-function(hills=hills) {
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

# head hills
head.hillsfile<-function(hills=hills, n=10) {
  return(head(hills$hillsfile, n=n))
}

# tail hills
tail.hillsfile<-function(hills=hills, n=10) {
  return(tail(hills$hillsfile, n=n))
}

# sum HILLS from Plumed
`+.hillsfile`<-function(hills1, hills2) {
  if(ncol(hills1$hillsfile)!=ncol(hills2$hillsfile)) {
    stop("you can sum only hills of same dimension")
  }
  if(hills1$per[1]!=hills2$per[1]) {
    stop("you can sum only hills of same periodicity")
  }
  if(ncol(hills1$hillsfile)==7 || ncol(hills1$hillsfile)==8) {
    if(hills1$per[2]!=hills2$per[2]) {
      stop("you can sum only hills of same periodicity")
    }
  }
  hills<-list(hillsfile=rbind(hills1$hillsfile, hills2$hillsfile), size=dim(rbind(hills1$hillsfile, hills2$hillsfile)),
              filename=hills1$filename, per=hills1$per, pcv1=hills1$pcv1, pcv2=hills1$pcv2)
  class(hills) <- "hillsfile"
  return(hills)
}

# plot hillsfile
plot.hillsfile<-function(hills=hills, ignoretime=FALSE,
                         xlab=NULL, ylab=NULL,
                         xlim=NULL, ylim=NULL,
                         main=NULL, sub=NULL,
                         pch=1, col="black", bg="red", cex=1,
                         asp=NULL, lwd=1, axes=TRUE) {
  xlims<-NULL
  ylims<-NULL
  if(!is.null(xlim)) {xlims<-xlim}
  if((hills$per[1]==T)&is.null(xlim)) {xlims<-hills$pcv1}
  if(!is.null(ylim)) {ylims<-ylim}
  if((hills$per[2]==T)&is.null(ylim)) {ylims<-hills$pcv2}
  if(hills$size[2]==5) {
    if(is.null(xlab)) xlab="time"
    if(is.null(ylab)) ylab="CV"
    if(ignoretime) {
      plot(seq(from=hills$hillsfile[1,1],by=hills$hillsfile[1,1],length.out=nrow(hills$hillsfile)),
           hills$hillsfile[,2], type="l",
           xlab=xlab, ylab=ylab,
           main=main, sub=sub,
           xlim=xlims, ylim=ylims,
           col=col, cex=cex, lwd=lwd,
           asp=asp, axes=axes)
    } else {
      plot(hills$hillsfile[,1], hills$hillsfile[,2], type="l",
           xlab=xlab, ylab=ylab,
           main=main, sub=sub,
           xlim=xlims, ylim=ylims,
           col=col, cex=cex, lwd=lwd,
           asp=asp, axes=axes)
    }
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

# points hillsfile
points.hillsfile<-function(hills=hills, ignoretime=FALSE,
                           pch=1, col="black", bg="red", cex=1,
                           asp=NULL, lwd=1, axes=TRUE) {
  if(hills$size[2]==5) {
    if(ignoretime) {
      points(seq(from=hills$hillsfile[1,1],by=hills$hillsfile[1,1],length.out=nrow(hills$hillsfile)),
             hills$hillsfile[,2],
             col=col, cex=cex, lwd=lwd)
    } else {
      points(hills$hillsfile[,1], hills$hillsfile[,2],
             col=col, cex=cex, lwd=lwd)
    }
  }
  if(hills$size[2]==7) {
    points(hills$hillsfile[,2], hills$hillsfile[,3],
           pch=pch, col=col, bg=bg, cex=cex, lwd=lwd)
  }
}

# lines hillsfile
lines.hillsfile<-function(hills=hills, ignoretime=FALSE,
                          lwd=1, col="black") {
  if(hills$size[2]==5) {
    if(ignoretime) {
      lines(seq(from=hills$hillsfile[1,1],by=hills$hillsfile[1,1],length.out=nrow(hills$hillsfile)),
            hills$hillsfile[,2],
            col=col, lwd=lwd)
    } else {
      plot(hills$hillsfile[,1], hills$hillsfile[,2],
           col=col, lwd=lwd)
    }
  }
  if(hills$size[2]==7) {
    lines(hills$hillsfile[,2], hills$hillsfile[,3],
          col=col, lwd=lwd)
  }
}

# plot heights
plotheights<-function(hills=hills, ignoretime=FALSE, xlab=NULL, ylab=NULL,
                      xlim=NULL, ylim=NULL, zlim=NULL,
                      main=NULL, sub=NULL,
                      pch=1, col="black", bg="red", cex=1,
                      asp=NULL, lwd=1, axes=TRUE) {
  if(class(hills)=="hillsfile") {
    if(is.null(xlab)) xlab="time"
    if(is.null(ylab)) ylab="hill height"
    if(hills$size[2]==5) {
      if(ignoretime) {
        plot(seq(from=hills$hillsfile[1,1],by=hills$hillsfile[1,1],length.out=nrow(hills$hillsfile)),
             hills$hillsfile[,4], type="l",
             xlab=xlab, ylab=ylab,
             main=main, sub=sub,
             col=col, cex=cex, lwd=lwd,
             asp=asp, axes=axes)
      } else {
        plot(hills$hillsfile[,1], hills$hillsfile[,4], type="l",
             xlab=xlab, ylab=ylab,
             main=main, sub=sub,
             col=col, cex=cex, lwd=lwd,
             asp=asp, axes=axes)
      }
    }
    if(hills$size[2]==7) {
      if(ignoretime) {
        plot(seq(from=hills$hillsfile[1,1],by=hills$hillsfile[1,1],length.out=nrow(hills$hillsfile)),
             hills$hillsfile[,6], type="l",
             xlab=xlab, ylab=ylab,
             main=main, sub=sub,
             col=col, cex=cex, lwd=lwd,
             asp=asp, axes=axes)
      } else {
        plot(hills$hillsfile[,1], hills$hillsfile[,6], type="l",
             xlab=xlab, ylab=ylab,
             main=main, sub=sub,
             col=col, cex=cex, lwd=lwd,
             asp=asp, axes=axes)
      }
    }
  } else {
    stop("function plotheights requires object hillsfile as an input")
  }
}

# red FES from MetadynView
read.fes<-function(filename=filename, dimension=2, per=c(TRUE, TRUE), pcv1=c(-pi,pi), pcv2=c(-pi,pi)) {
  ifile<-read.table(filename)
  rows<-sqrt(nrow(ifile))
  fes<-matrix(ifile[,3], nrow=rows)
  fes<-max(fes)-fes
  x<-min(ifile[,1])+(max(ifile[,1])-min(ifile[,1]))*0:(rows-1)/(rows-1)
  y<-min(ifile[,2])+(max(ifile[,2])-min(ifile[,2]))*0:(rows-1)/(rows-1)
  cfes<-list(fes=fes, rows=rows, dimension=dimension, per=per, x=x, y=y, pcv1=pcv1, pcv2=pcv2)
  class(cfes) <- "fes"
  return(cfes)
}

# calculate fes by bias sum algorithm
fes<-function(hills=hills, tmin=0, tmax=NULL, xlim=NULL, ylim=NULL, npoints=256) {
  if(!is.null(tmax)) {
    if(hills$size[1]<tmax) {
      cat("You requested more hills by tmax than available, using all hills\n")
      tmax<-hills$size[1]
    }
  }
  if(is.null(tmax)) {
    tmax<-hills$size[1]
  }
  if(tmin>=tmax) {
    stop("tmax must be higher than tmin")
  }
  sourceCpp("../src/mm.cpp")
  if(hills$size[2]==7) {
    if(max(hills$hillsfile[,4])/min(hills$hillsfile[,4])>1.00000000001) {
      stop("Bias Sum algorithm works only with hills of the same sizes")
    }
    if(max(hills$hillsfile[,5])/min(hills$hillsfile[,5])>1.00000000001) {
      stop("Bias Sum algorithm works only with hills of the same sizes")
    }
    minCV1 <- min(hills$hillsfile[,2])
    maxCV1 <- max(hills$hillsfile[,2])
    minCV2 <- min(hills$hillsfile[,3])
    maxCV2 <- max(hills$hillsfile[,3])
    xlims<-c(minCV1-0.05*(maxCV1-minCV1), maxCV1+0.05*(maxCV1-minCV1))
    ylims<-c(minCV2-0.05*(maxCV2-minCV2), maxCV2+0.05*(maxCV2-minCV2))
    if(!is.null(xlim)) {xlims<-xlim}
    if((hills$per[1]==T)&is.null(xlim)) {xlims<-hills$pcv1}
    if(!is.null(ylim)) {ylims<-ylim}
    if((hills$per[2]==T)&is.null(ylim)) {ylims<-hills$pcv2}
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
    cfes<-list(fes=fesm, hills=hills$hillsfile, rows=npoints, dimension=2, per=hills$per, x=x, y=y, pcv1=hills$pcv1, pcv2=hills$pcv2)
    class(cfes) <- "fes"
  }
  if(hills$size[2]==5) {
    if(max(hills$hillsfile[,3])/min(hills$hillsfile[,3])>1.00000000001) {
      stop("Bias Sum algorithm works only with hills of the same sizes")
    }
    minCV1 <- min(hills$hillsfile[,2])
    maxCV1 <- max(hills$hillsfile[,2])
    xlims<-c(minCV1-0.05*(maxCV1-minCV1), maxCV1+0.05*(maxCV1-minCV1))
    if(!is.null(xlim)) {xlims<-xlim}
    if((hills$per[1]==T)&is.null(xlim)) {xlims<-hills$pcv1}
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
    cfes<-list(fes=fesm, hills=hills$hillsfile, rows=npoints, dimension=1, per=hills$per, x=x, pcv1=hills$pcv1, pcv2=hills$pcv2)
    class(cfes) <- "fes"
  }
  return(cfes)
}

# calculate fes conventionally (slow)
fes2<-function(hills=hills, tmin=0, tmax=NULL, xlim=NULL, ylim=NULL, npoints=256) {
  if(!is.null(tmax)) {
    if(hills$size[1]<tmax) {
      cat("You requested more hills by tmax than available, using all hills\n")
      tmax<-hills$size[1]
    }
  }
  if(is.null(tmax)) {
    tmax<-hills$size[1]
  }
  if(tmin>=tmax) {
    stop("tmax must be higher than tmin")
  }
  sourceCpp("../src/mm.cpp")
  if(hills$size[2]==7) {
    minCV1 <- min(hills$hillsfile[,2])
    maxCV1 <- max(hills$hillsfile[,2])
    minCV2 <- min(hills$hillsfile[,3])
    maxCV2 <- max(hills$hillsfile[,3])
    xlims<-c(minCV1-0.05*(maxCV1-minCV1), maxCV1+0.05*(maxCV1-minCV1))
    ylims<-c(minCV2-0.05*(maxCV2-minCV2), maxCV2+0.05*(maxCV2-minCV2))
    if(!is.null(xlim)) {xlims<-xlim}
    if((hills$per[1]==T)&is.null(xlim)) {xlims<-hills$pcv1}
    if(!is.null(ylim)) {ylims<-ylim}
    if((hills$per[2]==T)&is.null(ylim)) {ylims<-hills$pcv2}
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
    cfes<-list(fes=fesm, hills=hills$hillsfile, rows=npoints, dimension=2, per=hills$per, x=x, y=y, pcv1=hills$pcv1, pcv2=hills$pcv2)
    class(cfes) <- "fes"
  }
  if(hills$size[2]==5) {
    minCV1 <- min(hills$hillsfile[,2])
    maxCV1 <- max(hills$hillsfile[,2])
    xlims<-c(minCV1-0.05*(maxCV1-minCV1), maxCV1+0.05*(maxCV1-minCV1))
    if(!is.null(xlim)) {xlims<-xlim}
    if((hills$per[1]==T)&is.null(xlim)) {xlims<-hills$pcv1}
    x<-0:(npoints-1)*(xlims[2]-xlims[1])/(npoints-1)+xlims[1]
    if(hills$per[1]==F) {
      fesm<-hills1d2(npoints*(hills$hillsfile[,2]-xlims[1])/(xlims[2]-xlims[1]),
                     npoints*hills$hillsfile[,3]/(xlims[2]-xlims[1]),
                     hills$hillsfile[,4],npoints,tmin,tmax)
    }
    if(hills$per[1]==T) {
      fesm<-hills1d2p(npoints*(hills$hillsfile[,2]-xlims[1])/(xlims[2]-xlims[1]),
                      npoints*hills$hillsfile[,3]/(xlims[2]-xlims[1]),
                      hills$hillsfile[,4],npoints,tmin,tmax)
    }
    cfes<-list(fes=fesm, hills=hills$hillsfile, rows=npoints, dimension=1, per=hills$per, x=x, pcv1=hills$pcv1, pcv2=hills$pcv2)
    class(cfes) <- "fes"
  }
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
      cfes<-list(fes=fes1$fes+fes2$fes, hills=rbind(fes1$hillsfile, fes2$hillsfile), rows=fes1$rows, dimension=fes1$dimension, per=fes1$per, x=fes1$x, pcv1=fes1$pcv1, pcv2=fes1$pcv2)
    }
    if(fes1$dimension==2) {
      cfes<-list(fes=fes1$fes+fes2$fes, hills=rbind(fes1$hillsfile, fes2$hillsfile), rows=fes1$rows, dimension=fes1$dimension, per=fes1$per, x=fes1$x, y=fes1$y, pcv1=fes1$pcv1, pcv2=fes1$pcv2)
    }
  } else if(class(fes1)=="fes") {
    if(fes1$dimension==1) {
      cfes<-list(fes=fes1$fes+fes2, hills=fes1$hillsfile, rows=fes1$rows, dimension=fes1$dimension, per=fes1$per, x=fes1$x, pcv1=fes1$pcv1, pcv2=fes1$pcv2)
    }
    if(fes1$dimension==2) {
      cfes<-list(fes=fes1$fes+fes2, hills=fes1$hillsfile, rows=fes1$rows, dimension=fes1$dimension, per=fes1$per, x=fes1$x, y=fes1$y, pcv1=fes1$pcv1, pcv2=fes1$pcv2)
    }
  } else if(class(fes2)=="fes") {
    if(fes2$dimension==1) {
      cfes<-list(fes=fes1+fes2$fes, hills=fes2$hillsfile, rows=fes2$rows, dimension=fes2$dimension, per=fes2$per, x=fes2$x, pcv1=fes2$pcv1, pcv2=fes2$pcv2)
    }
    if(fes2$dimension==2) {
      cfes<-list(fes=fes1+fes2$fes, hills=rbind(fes1$hillsfile,fes2$hillsfile), rows=fes2$rows, dimension=fes2$dimension, per=fes2$per, x=fes2$x, y=fes2$y, pcv1=fes2$pcv1, pcv2=fes2$pcv2)
    }
  }
  class(cfes) <- "fes"
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
    cat("WARNING: FES obtained by subtraction of two FESes\n")
    cat(" will inherit hills only from the first FES\n")
    if(fes1$dimension==1) {
      cfes<-list(fes=fes1$fes-fes2$fes, hills=fes1$hillsfile, rows=fes1$rows, dimension=fes1$dimension, per=fes1$per, x=fes1$x, pcv1=fes1$pcv1, pcv2=fes1$pcv2)
    }
    if(fes1$dimension==2) {
      cfes<-list(fes=fes1$fes-fes2$fes, hills=fes1$hillsfile, rows=fes1$rows, dimension=fes1$dimension, per=fes1$per, x=fes1$x, y=fes1$y, pcv1=fes1$pcv1, pcv2=fes1$pcv2)
    }
  } else if(class(fes1)=="fes") {
    if(fes1$dimension==1) {
      cfes<-list(fes=fes1$fes-fes2, hills=fes1$hillsfile, rows=fes1$rows, dimension=fes1$dimension, per=fes1$per, x=fes1$x, pcv1=fes1$pcv1, pcv2=fes1$pcv2)
    }
    if(fes1$dimension==2) {
      cfes<-list(fes=fes1$fes-fes2, hills=fes1$hillsfile, rows=fes1$rows, dimension=fes1$dimension, per=fes1$per, x=fes1$x, y=fes1$y, pcv1=fes1$pcv1, pcv2=fes1$pcv2)
    }
  } else if(class(fes2)=="fes") {
    if(fes2$dimension==1) {
      cfes<-list(fes=fes1-fes2$fes, hills=fes2$hillsfile, rows=fes2$rows, dimension=fes2$dimension, per=fes2$per, x=fes2$x, pcv1=fes2$pcv1, pcv2=fes2$pcv2)
    }
    if(fes2$dimension==2) {
      cfes<-list(fes=fes1-fes2$fes, hills=fes2$hillsfile, rows=fes2$rows, dimension=fes2$dimension, per=fes2$per, x=fes2$x, y=fes2$y, pcv1=fes2$pcv1, pcv2=fes2$pcv2)
    }
  }
  class(cfes) <- "fes"
  return(cfes)
}

# multiply a fes
`*.fes`<-function(fes1, fes2) {
  if((class(fes1)=="fes")&(class(fes2)=="fes")) {
    stop("you cannot multiply fes by fes")
  } else if(class(fes1)=="fes") {
    if(fes1$dimension==1) {
      cfes<-list(fes=fes1$fes*fes2, hills=fes1$hillsfile, rows=fes1$rows, dimension=fes1$dimension, per=fes1$per, x=fes1$x, pcv1=fes1$pcv1, pcv2=fes1$pcv2)
    }
    if(fes1$dimension==2) {
      cfes<-list(fes=fes1$fes*fes2, hills=fes1$hillsfile, rows=fes1$rows, dimension=fes1$dimension, per=fes1$per, x=fes1$x, y=fes1$y, pcv1=fes1$pcv1, pcv2=fes1$pcv2)
    }
  } else if(class(fes2)=="fes") {
    if(fes2$dimension==1) {
      cfes<-list(fes=fes1*fes2$fes, hills=fes2$hillsfile, rows=fes2$rows, dimension=fes2$dimension, per=fes2$per, x=fes2$x, pcv1=fes2$pcv1, pcv2=fes2$pcv2)
    }
    if(fes2$dimension==2) {
      cfes<-list(fes=fes1*fes2$fes, hills=fes2$hillsfile, rows=fes2$rows, dimension=fes2$dimension, per=fes2$per, x=fes2$x, y=fes2$y, pcv1=fes2$pcv1, pcv2=fes2$pcv2)
    }
  }
  cat("WARNING: multiplication of FES will multiply\n")
  cat(" the FES but not hill heights\n")
  class(cfes) <- "fes"
  return(cfes)
}

# divide a fes
`/.fes`<-function(fes1, coef) {
  if((class(fes1)=="fes")&(class(coef)=="fes")) {
    stop("you cannot divide fes by fes")
  } else if(class(fes1)=="fes") {
    if(fes1$dimension==1) {
      cfes<-list(fes=fes1$fes/coef, hills=fes1$hillsfile, rows=fes1$rows, dimension=fes1$dimension, per=fes1$per, x=fes1$x, pcv1=fes1$pcv1, pcv2=fes1$pcv2)
    }
    if(fes1$dimension==2) {
      cfes<-list(fes=fes1$fes/coef, hills=fes1$hillsfile, rows=fes1$rows, dimension=fes1$dimension, per=fes1$per, x=fes1$x, y=fes1$y, pcv1=fes1$pcv1, pcv2=fes1$pcv2)
    }
  } else if(class(coef)=="fes") {
    stop("you cannot divide something by fes")
  }
  cat("WARNING: division of FES will divide\n")
  cat(" the FES but not hill heights\n")
  class(cfes) <- "fes"
  return(cfes)
}

# min of fes
min.fes<-function(inputfes=inputfes, na.rm=NULL) {
  return(min(inputfes$fes, na.rm=na.rm))
}

# max of fes
max.fes<-function(inputfes=inputfes, na.rm=NULL) {
  return(max(inputfes$fes, na.rm=na.rm))
}

# mean of fes
mean.fes<-function(inputfes=inputfes, na.rm=NULL) {
  return(mean(inputfes$fes, na.rm=na.rm))
}

# print FES
print.fes<-function(inputfes=inputfes) {
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
summary.fes<-function(inputfes=inputfes) {
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
plot.fes<-function(inputfes=inputfes, plottype="both",
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

# points FES
points.fes<-function(inputfes=inputfes, x=NULL,
                     pch=1, col="black", bg="red", cex=1) {
  fes<-inputfes$fes
  if(inputfes$dimension==1) {
    if(is.null(x)) x<-inputfes$x
    points(x, fes,
           pch=pch, col=col, bg=bg, cex=cex)
  } else {
    cat("points available only for 1D free energy surfaces\n")
  }
}

# lines FES
lines.fes<-function(inputfes=inputfes, x=NULL,
                    lwd=1, col="black") {
  fes<-inputfes$fes
  if(inputfes$dimension==1) {
    if(is.null(x)) x<-inputfes$x
    lines(x, fes, lwd=lwd, col=col)
  } else {
    cat("points available only for 1D free energy surfaces\n")
  }
}

# find minima of a FES
fesminima<-function(inputfes=inputfes, nbins=8) {
  fes<-inputfes$fes
  rows<-inputfes$rows
  rb <- rows/nbins
  if(rb<2) {
    stop("nbins too high, try to reduce it")
  }
  if(rows%%nbins>0) {
    stop("number of rows in FES must be integer multiple of nbins")
  }
  per<-inputfes$per
  if(inputfes$dimension==2) {
    minx<-c()
    miny<-c()
    for(i in 0:(nbins-1)) {
      ni<-i*rb+0:(rb+1)
      if(per[1]) {
        ni[ni==0]<-rows
        ni[ni==(rows+1)]<-1
      } else {
        ni<-ni[ni!=0]
        ni<-ni[ni!=(rows+1)]
      }
      for(j in 0:(nbins-1)) {
        nj<-j*rb+0:(rb+1)
        if(per[2]) {
          nj[nj==0]<-rows
          nj[nj==(rows+1)]<-1
        } else {
          nj<-nj[nj!=0]
          nj<-nj[nj!=(rows+1)]
        }
        binmin<-which(fes[ni,nj]==min(fes[ni,nj]), arr.ind = TRUE)
        if(binmin[1]!=1 && binmin[2]!=1 && binmin[1]!=length(ni) && binmin[2]!=length(nj)) {
          minx<-c(minx,i*rb+binmin[1]-1)
          miny<-c(miny,j*rb+binmin[2]-1)
        }
      }
    }
    myLETTERS <- c(LETTERS, paste("A", LETTERS, sep=""), paste("B", LETTERS, sep=""))[1:length(minx)]
    minima<-data.frame(myLETTERS, minx, miny, inputfes$x[minx], inputfes$y[miny], fes[cbind(minx,miny)])
    names(minima) <- c("letter", "CV1bin", "CV2bin", "CV1", "CV2", "free_energy")
    minima <- minima[order(minima[,6]),]
    rownames(minima) <- seq(length=nrow(minima))
    minima[,1]<-myLETTERS
    minima<-list(minima=minima, hills=inputfes$hills, fes=fes, rows=rows, dimension=inputfes$dimension, per=per, x=inputfes$x, y=inputfes$y, pcv1=inputfes$pcv1, pcv2=inputfes$pcv2)
    class(minima) <- "minima"
  }
  if(inputfes$dimension==1) {
    minx<-c()
    for(i in 0:(nbins-1)) {
      ni<-i*rb+0:(rb+1)
      if(per[1]) {
        ni[ni==0]<-rows
        ni[ni==(rows+1)]<-1
      } else {
        ni<-ni[ni!=0]
        ni<-ni[ni!=(rows+1)]
      }
      binmin<-which(fes[ni]==min(fes[ni]), arr.ind = TRUE)
      if(binmin[1]!=1 && binmin[1]!=length(ni)) {
        minx<-c(minx,i*rb+binmin[1]-1)
      }
    }
    myLETTERS <- c(LETTERS, paste("A", LETTERS, sep=""), paste("B", LETTERS, sep=""))[1:length(minx)]
    minima<-data.frame(myLETTERS, minx, inputfes$x[minx], fes[minx])
    names(minima) <- c("letter", "CV1bin", "CV1", "free_energy")
    minima <- minima[order(minima[,4]),]
    rownames(minima) <- seq(length=nrow(minima))
    minima[,1]<-myLETTERS
    minima<-list(minima=minima, hills=inputfes$hills, fes=fes, rows=rows, dimension=inputfes$dimension, per=per, x=inputfes$x, pcv1=inputfes$pcv1, pcv2=inputfes$pcv2)
    class(minima) <- "minima"
  }
  return(minima)
}

# create empty minima
emptyminima<-function(inputfes=inputfes) {
  fes<-inputfes$fes
  rows<-inputfes$rows
  per<-inputfes$per
  if(inputfes$dimension==2) {
    minima<-data.frame(c("A"), c(0), c(0), c(0), c(0), c(0))
    minima<-minima[-1,]
    names(minima) <- c("letter", "CV1bin", "CV2bin", "CV1", "CV2", "free_energy")
    minima<-list(minima=minima, hills=inputfes$hills, fes=fes, rows=rows, dimension=inputfes$dimension, per=per, x=inputfes$x, y=inputfes$y, pcv1=inputfes$pcv1, pcv2=inputfes$pcv2)
    class(minima) <- "minima"
  }
  if(inputfes$dimension==1) {
    minima<-data.frame(c("A"), c(0), c(0), c(0))
    minima<-minima[-1,]
    names(minima) <- c("letter", "CV1bin", "CV1", "free_energy")
    minima<-list(minima=minima, hills=inputfes$hills, fes=fes, rows=rows, dimension=inputfes$dimension, per=per, x=inputfes$x, pcv1=inputfes$pcv1, pcv2=inputfes$pcv2)
    class(minima) <- "minima"
  }
  return(minima)
}

# create one minima
oneminimum<-function(inputfes=inputfes, cv1=cv1, cv2=cv2) {
  fes<-inputfes$fes
  rows<-inputfes$rows
  per<-inputfes$per
  if(inputfes$dimension==2) {
    icv1<-as.integer(rows*(cv1-min(inputfes$x))/(max(inputfes$x)-min(inputfes$x)))+1
    if(icv1<0)    stop("out of range")
    if(icv1>rows) stop("out of range")
    icv2<-as.integer(rows*(cv2-min(inputfes$y))/(max(inputfes$x)-min(inputfes$x)))+1
    if(icv2<0)    stop("out of range")
    if(icv2>rows) stop("out of range")
    minima<-data.frame(c("A"), c(icv1), c(icv2), c(cv1), c(cv2), c(fes[icv1,icv2]))
    names(minima) <- c("letter", "CV1bin", "CV2bin", "CV1", "CV2", "free_energy")
    minima<-list(minima=minima, hills=inputfes$hills, fes=fes, rows=rows, dimension=inputfes$dimension, per=per, x=inputfes$x, y=inputfes$y, pcv1=inputfes$pcv1, pcv2=inputfes$pcv2)
    class(minima) <- "minima"
  }
  if(inputfes$dimension==1) {
    icv1<-as.integer(rows*(cv1-min(inputfes$x))/(max(inputfes$x)-min(inputfes$x)))+1
    if(icv1<0)    stop("out of range")
    if(icv1>rows) stop("out of range")
    minima<-data.frame(c("A"), c(icv1), c(cv1), c(fes[icv1,icv2]))
    names(minima) <- c("letter", "CV1bin", "CV1", "free_energy")
    minima<-list(minima=minima, hills=inputfes$hills, fes=fes, rows=rows, dimension=inputfes$dimension, per=per, x=inputfes$x, pcv1=inputfes$pcv1, pcv2=inputfes$pcv2)
    class(minima) <- "minima"
  }
  return(minima)
}

# add minima
`+.minima`<-function(min1, min2) {
  if(class(min1)!="minima") {
    stop("you can sum only two minima objects")
  }
  if(class(min2)!="minima") {
    stop("you can sum only two minima objects")
  }
  if(sum(min1$fes)!=sum(min2$fes)) {
    stop("you can sum only minima objects with same FESes")
  }
  myLETTERS <- c(LETTERS, paste("A", LETTERS, sep=""), paste("B", LETTERS, sep=""))[1:(nrow(min1$minima)+nrow(min2$minima))]
  minima1<-min1$minima
  minima2<-min2$minima
  minima<-rbind(minima1, minima2)
  if(inputfes$dimension==2) {
    names(minima) <- c("letter", "CV1bin", "CV2bin", "CV1", "CV2", "free_energy")
    minima <- minima[order(minima[,6]),]
    rownames(minima) <- seq(length=nrow(minima))
    minima[,1]<-myLETTERS
    minima<-list(minima=minima, hills=min1$hills, fes=min1$fes, rows=min1$rows, dimension=min1$dimension, per=min1$per, x=min1$x, y=min1$y, pcv1=min1$pcv1, pcv2=min1$pcv2)
    class(minima) <- "minima"
  }
  if(inputfes$dimension==1) {
    names(minima) <- c("letter", "CV1bin", "CV1", "free_energy")
    minima <- minima[order(minima[,4]),]
    rownames(minima) <- seq(length=nrow(minima))
    minima[,1]<-myLETTERS
    minima<-list(minima=minima, hills=min1$hillsfile, fes=min1$fes, rows=min1$rows, dimension=min1$dimension, per=min1$per, x=min1$x, pcv1=min1$pcv1, pcv2=min1$pcv2)
    class(minima) <- "minima"
  }
  return(minima)
}

# print minima of a FES
print.minima<-function(minims) {
  cat("$minima\n\n")
  print(minims$minima)
}

# print a summary of minima of a FES
summary.minima<-function(minims=minims, temp=300, eunit="kJ/mol") {
  toprint <- minims$minima
  tind = 6
  if(minims$dimension==1) {
    tind = 4
  }
  print(tind)
  if(eunit=="kJ/mol") {
    toprint<-cbind(toprint, exp(-1000*toprint[,tind]/8.314/temp))
  }
  if(eunit=="J/mol") {
    toprint<-cbind(toprint, exp(-toprint[,tind]/8.314/temp))
  }
  if(eunit=="kcal/mol") {
    toprint<-cbind(toprint, exp(-1000*toprint[,tind]/8.314/temp/4.184))
  }
  if(eunit=="cal/mol") {
    toprint<-cbind(toprint, exp(-toprint[,tind]/8.314/temp/4.184))
  }
  sumpop<-sum(toprint[,tind+1])
  toprint<-cbind(toprint, 100*toprint[,tind+1]/sumpop)
  names(toprint)[tind+1]<-"relative_pop"
  names(toprint)[tind+2]<-"pop"
  print(toprint)
}

# plot minima
plot.minima <- function(minims=minims, plottype="both",
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
  minlabs<-minims$minima[,1]
  if(minims$dimension==1) {
    if(is.null(x)) x<-minims$x
    if(is.null(xlab)) xlab="CV"
    if(is.null(ylab)) ylab="free energy"
    if(is.null(xlim)) xlim<-c(min(x),max(x))
    if(is.null(ylim)) {
      ylim<-range(pretty(range(fes)))
    }
    minpoints<-minims$minima[,3:4]
    minpoints[,2]<-minpoints[,2]+0.05*(ylim[2]-ylim[1])
    plot(x, fes, type="l", lwd=lwd,
        col=col, xlim=xlim, ylim=ylim,
        xlab=xlab, ylab=ylab, axes=axes,
        main=main, sub=sub)
    text(minpoints, labels=minlabs, col=textcol, xlim=xlim, ylim=ylim, cex=cex)
  } else {
    minpoints<-minims$minima[,4:5]
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
      text(minpoints, labels=minlabs, col=textcol, xlim=xlim, ylim=ylim, cex=cex)
    }
    if(plottype=="contour") {
      contour(x, y, fes, zlim=zlim,
              nlevels=nlevels, levels=levels,
              labels=labels, labcex=labcex, drawlabels=drawlabels,
              method=method, col=contcol, lty=lty, lwd=lwd,
              main=main, sub=sub)
      text(minpoints, labels=minlabs, col=textcol, xlim=xlim, ylim=ylim, cex=cex)
    }
    if(plottype=="both") {
      contour(x, y, fes, zlim=zlim,
              nlevels=nlevels, levels=levels,
              labels=labels, labcex=labcex, drawlabels=drawlabels,
              method=method, col=contcol, lty=lty, lwd=lwd, add=T)
    }
  }
}

# Calculate free energy profiles
feprof <- function(minims=minims, tmin=0, tmax=NULL) {
  fes<-minims$fes
  rows<-minims$rows
  mins<-minims$minima
  hills<-minims$hills
  if(is.null(tmax)) {
    tmax<-nrow(hills)
  }
  if(tmax>nrow(hills)) {
    tmax<-nrow(hills)
    cat("You requested more hills by tmax than available, using all hills\n")
  }
  if(tmin>=tmax) {
    stop("tmax must be higher than tmin")
  }
  tt <- tmin:tmax
  mms <- data.frame(tt)
  if(minims$dimension==1) {
    for(i in 1:nrow(mins)) {
      if(minims$per[1]==T) {
        mm<-fe1dp(hills[,2], hills[,3], hills[,4], mins[i,3], minims$pcv1[2]-minims$pcv1[1], tmin, tmax)
      } else {
        mm<-fe1d(hills[,2], hills[,3], hills[,4], mins[i,3], tmin, tmax)
      }
      mms<-cbind(mms,mm)
    }
  }
  if(minims$dimension==2) {
    for(i in 1:nrow(mins)) {
      if(minims$per[1]==T && minims$per[2]==T) {
        mm<-fe2dp12(hills[,2], hills[,3], hills[,4], hills[,5], hills[,6], mins[i,4], mins[i,5], minims$pcv1[2]-minims$pcv1[1], minims$pcv2[2]-minims$pcv2[1], tmin, tmax)
      }
      if(minims$per[1]==T && minims$per[2]==F) {
        mm<-fe2dp1(hills[,2], hills[,3], hills[,4], hills[,5], hills[,6], mins[i,4], mins[i,5], minims$pcv1[2]-minims$pcv1[1], tmin, tmax)
      }
      if(minims$per[1]==F && minims$per[2]==T) {
        mm<-fe2dp2(hills[,2], hills[,3], hills[,4], hills[,5], hills[,6], mins[i,4], mins[i,5], minims$pcv2[2]-minims$pcv2[1], tmin, tmax)
      }
      if(minims$per[1]==F && minims$per[2]==F) {
        mm<-fe2d(hills[,2], hills[,3], hills[,4], hills[,5], hills[,6], mins[i,4], mins[i,5], tmin, tmax)
      }
      mms<-cbind(mms,mm)
    }
  }
  profs<-list(mms=mms, mins=mins, fes=fes, rows=rows, dimension=minims$dimension, per=minims$per, pcv1=minims$pcv1, pcv2=minims$pcv2)
  class(profs) <- "profiles"
  return(profs)
}

print.profiles <- function(profs=profs) {
  cat(nrow(profs$mms))
  if(profs$dimension==1) {
    cat(" 1D minima\n")
  } else {
    cat(" 2D minima\n")
  }
}

summary.profiles <- function(profs=profs) {
  if(profs$dimension==1) {
    outprofile <- profs$mins
    mms<-profs$mms[,2:ncol(profs$mms)]-profs$mms[,2]
    outprofile <- cbind(outprofile,apply(mms,2,min))
    outprofile <- cbind(outprofile,apply(mms,2,max))
    outprofile <- cbind(outprofile,t(mms[nrow(mms),]))
    names(outprofile)[5:7]<-c("min diff", "max diff", "tail")
    print(outprofile)
  } else {
    outprofile <- profs$mins
    mms<-profs$mms[,2:ncol(profs$mms)]-profs$mms[,2]
    outprofile <- cbind(outprofile,apply(mms,2,min))
    outprofile <- cbind(outprofile,apply(mms,2,max))
    outprofile <- cbind(outprofile,t(mms[nrow(mms),]))
    names(outprofile)[7:9]<-c("min diff", "max diff", "tail")
    print(outprofile)
  }
}

plot.profiles <- function(profs=profs, which=NULL,
                          ignoretime=FALSE,
                          xlim=NULL, ylim=NULL,
                          main=NULL, sub=NULL,
                          xlab=NULL, ylab=NULL,
                          col=NULL, asp=NULL, lwd=1, axes=T) {
  if(is.null(which)) which<-1:(ncol(profs$mms)-1)
  if(is.null(xlab)) xlab<-"Index"
  if(is.null(ylab)) ylab<-"Free Energy Difference (kJ/mol)"
  if(is.null(col)) col<-rainbow(ceiling(1.35*length(which)))[length(which):1]
  col<-rep(col, (ncol(profs$mms)-1))
  mms<-profs$mms[,2:ncol(profs$mms)]-profs$mms[,2]
  if(is.null(xlim)) xlim<-c(min(profs$mms[,1]),max(profs$mms[,1]))
  if(is.null(ylim)) {
    yliml <- min(mms[,1:ncol(mms)])
    ylimu <- max(mms[,1:ncol(mms)])
    ylim<-c(yliml-0.05*(ylimu-yliml), ylimu+0.05*(ylimu-yliml))
  }
  if(ignoretime) {
    plot(seq(from=profs$mms[1,1], by=profs$mms[1,1], length.out=nrow(profs$mms)),
         mms[,1], type="l",
         xlab=xlab, ylab=ylab,
         main=main, sub=sub,
         xlim=xlim, ylim=ylim,
         lwd=lwd, asp=asp, axes=axes)
    for(i in 1:length(which)) {
        lines(profs$mms[,1], mms[,which[i]],
              lwd=lwd, col=col[i])
    }
  } else {
    plot(profs$mms[,1], mms[,1], type="l",
         xlab=xlab, ylab=ylab,
         main=main, sub=sub,
         xlim=xlim, ylim=ylim,
         lwd=lwd, asp=asp, axes=axes)
    for(i in 1:length(which)) {
        lines(profs$mms[,1], mms[,which[i]],
              lwd=lwd, col=col[i])
    }
  }
}

neb<-function(minims, min1, min2, nbins=20,
              nsteps=100, step=1.0, k=0.2) {
  fes<-minims$fes
  pcv1<-minims$pcv1
  pcv2<-minims$pcv2
  rows<-minims$rows
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
       f1 <- force1(fes, x, y)
       f1a <- align(f1, tau)
       f1p <- f1-f1a
       f2 <- force2(path[j,],path[j+1,],path[j+2,],k)
       f2a <- align(f2, tau)
       f2p <- f2-f2a
       fphi <- 0.5*(1.0+cos(cosphi(path[j,], path[j+1,], path[j+2,])*pi))
       newpath[j+1,1] <- newpath[j+1,1] + step*(f1p[1]+f2a[1]+fphi*f2p[1])
       newpath[j+1,2] <- newpath[j+1,2] + step*(f1p[2]+f2a[2]+fphi*f2p[2])
       newpath[j+1,3] <- fes[newpath[j+1,1],newpath[j+1,2]]
    }
    path <- newpath
  }
  path[,1]<-path[,1]*(pcv1[2]-pcv1[1])/rows+pcv1[1]
  path[,2]<-path[,2]*(pcv2[2]-pcv2[1])/rows+pcv2[1]
  #profs<-list(mms=mms, mins=mins, fes=fes, rows=rows, dimension=minims$dimension, per=minims$per, pcv1=minims$pcv1, pcv2=minims$pcv2)
  #class(profs) <- "profiles"
  #return(profs)
  #make object
  return(path)
}

# plot path
# print path
# summary path
# min path
# max path


