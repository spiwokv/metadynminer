#' Find free energy minima in the fes object
#'
#' `fesminima` finds free energy minima on 1D or 2D free energy surface.
#' The surface is divided by a 1D or 2D grid and minima are found for each
#' grid point. Next the program determines whether the minimum is a local
#' free energy minimum. Free energy minima are labeled constitutively by
#' capital letters.
#'
#' @param inputfes fes object
#' @param nbins number of bins for each CV (default 8)
#' @return minima object
#'
#' @export
#' @examples
#' tfes<-fes(acealanme)
#' minima<-fesminima(tfes)
#' minima
fesminima<-function(inputfes=inputfes, nbins=8) {
  fes<-inputfes$fes
  rows<-inputfes$rows
  rb <- rows/nbins
  if(rb<2) {
    stop("Error: nbins too high, try to reduce it")
  }
  if(rows%%nbins>0) {
    stop("Error: number of rows in FES must be integer multiple of nbins")
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

#' @export
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

#' Creates one ad hoc free energy minimum for a fes object
#'
#' `oneminimum` creates an ad hoc free energy minimum on free energy surface.
#' This can be used to calculate free energy surface evolution at arbitrary
#' point of free energy surface.
#'
#' @param inputfes fes object
#' @param cv1 the value of collective variable 1
#' @param cv2 the value of collective variable 2
#' @return minima object
#'
#' @export
#' @examples
#' tfes<-fes(acealanme1d)
#' minima<-fesminima(tfes)
#' minima<-minima+oneminimum(tfes, cv1=0, cv2=0)
#' minima
oneminimum<-function(inputfes=inputfes, cv1=cv1, cv2=cv2) {
  fes<-inputfes$fes
  rows<-inputfes$rows
  per<-inputfes$per
  if(inputfes$dimension==2) {
    icv1<-as.integer(rows*(cv1-min(inputfes$x))/(max(inputfes$x)-min(inputfes$x)))+1
    if(icv1<0)    stop("Error: Out of range")
    if(icv1>rows) stop("Error: Out of range")
    icv2<-as.integer(rows*(cv2-min(inputfes$y))/(max(inputfes$x)-min(inputfes$x)))+1
    if(icv2<0)    stop("Error: Out of range")
    if(icv2>rows) stop("Error: Out of range")
    minima<-data.frame(c("A"), c(icv1), c(icv2), c(cv1), c(cv2), c(fes[icv1,icv2]))
    names(minima) <- c("letter", "CV1bin", "CV2bin", "CV1", "CV2", "free_energy")
    minima<-list(minima=minima, hills=inputfes$hills, fes=fes, rows=rows, dimension=inputfes$dimension, per=per, x=inputfes$x, y=inputfes$y, pcv1=inputfes$pcv1, pcv2=inputfes$pcv2)
    class(minima) <- "minima"
  }
  if(inputfes$dimension==1) {
    icv1<-as.integer(rows*(cv1-min(inputfes$x))/(max(inputfes$x)-min(inputfes$x)))+1
    if(icv1<0)    stop("Error: Out of range")
    if(icv1>rows) stop("Error: Out of range")
    minima<-data.frame(c("A"), c(icv1), c(cv1), c(fes[icv1]))
    names(minima) <- c("letter", "CV1bin", "CV1", "free_energy")
    minima<-list(minima=minima, hills=inputfes$hills, fes=fes, rows=rows, dimension=inputfes$dimension, per=per, x=inputfes$x, pcv1=inputfes$pcv1, pcv2=inputfes$pcv2)
    class(minima) <- "minima"
  }
  return(minima)
}

#' @export
`+.minima`<-function(min1, min2) {
  if(class(min1)!="minima") {
    stop("Error: You can sum only two minima objects")
  }
  if(class(min2)!="minima") {
    stop("Error: You can sum only two minima objects")
  }
  if(sum(min1$fes)!=sum(min2$fes)) {
    stop("Error: You can sum only minima objects with same FESes")
  }
  myLETTERS <- c(LETTERS, paste("A", LETTERS, sep=""), paste("B", LETTERS, sep=""))[1:(nrow(min1$minima)+nrow(min2$minima))]
  minima1<-min1$minima
  minima2<-min2$minima
  minima<-rbind(minima1, minima2)
  if(min1$dimension==2) {
    names(minima) <- c("letter", "CV1bin", "CV2bin", "CV1", "CV2", "free_energy")
    minima <- minima[order(minima[,6]),]
    rownames(minima) <- seq(length=nrow(minima))
    minima[,1]<-myLETTERS
    minima<-list(minima=minima, hills=min1$hills, fes=min1$fes, rows=min1$rows, dimension=min1$dimension, per=min1$per, x=min1$x, y=min1$y, pcv1=min1$pcv1, pcv2=min1$pcv2)
    class(minima) <- "minima"
  }
  if(min1$dimension==1) {
    names(minima) <- c("letter", "CV1bin", "CV1", "free_energy")
    minima <- minima[order(minima[,4]),]
    rownames(minima) <- seq(length=nrow(minima))
    minima[,1]<-myLETTERS
    minima<-list(minima=minima, hills=min1$hillsfile, fes=min1$fes, rows=min1$rows, dimension=min1$dimension, per=min1$per, x=min1$x, pcv1=min1$pcv1, pcv2=min1$pcv2)
    class(minima) <- "minima"
  }
  return(minima)
}

#' Print minima object
#'
#' `print.minima` prints free energy minima (identifier, values of bins and collective variables and free energy).
#'
#' @param minims minima object
#' @param temp temperature in Kelvins
#' @param eunit energy units (kJ/mol or kcal/mol, kJ/mol is default)
#'
#' @export
#' @examples
#' tfes<-fes(acealanme)
#' minima<-fesminima(tfes)
#' minima
print.minima<-function(minims) {
  cat("$minima\n\n")
  print(minims$minima)
}

#' Print minima object summary
#'
#' `summary.minima` prints summary for free energy minima (identifier, values of bins and collective variables,
#' free energy and equilibrium populations).
#'
#' @param minims minima object
#' @param temp temperature in Kelvins
#' @param eunit energy units (kJ/mol or kcal/mol, kJ/mol is default)
#'
#' @export
#' @examples
#' tfes<-fes(acealanme)
#' minima<-fesminima(tfes)
#' summary(minima)
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
  if(eunit=="kcal/mol") {
    toprint<-cbind(toprint, exp(-1000*toprint[,tind]/8.314/temp/4.184))
  }
  sumpop<-sum(toprint[,tind+1])
  toprint<-cbind(toprint, 100*toprint[,tind+1]/sumpop)
  names(toprint)[tind+1]<-"relative_pop"
  names(toprint)[tind+2]<-"pop"
  print(toprint)
}

#' Plot minima object
#'
#' `plot.minima` plots free energy surface with minima. The free energy sufrace is ploted the same
#' way as by plot.fes with aditional minima labels.
#'
#' @param minims minima object
#' @param plottype specifies whether 2D free energy surface will be ploted as image, contours or both (default "both")
#' @param colscale specifies whether color scale will be ploted (default False)
#' @param colscalelab color scale label (default "free energy")
#' @inherit plot
#' @inherit image
#' @inherit contours
#'
#' @export
#' @examples
#' tfes<-fes(acealanme)
#' minima<-fesminima(tfes)
#' plot(minima)
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
                  axes=TRUE) {
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

#' Calculate free energy profile for minima object
#'
#' `feprof` calculates free energy profiles for free energy minima. It finds the global minimum at the
#' imax and calculates the evolution of free of local vs. global free energy minimum. The global
#' minimum is const (zero).
#'
#' @param minims minima object
#' @param imin index of a hill from which summation starts (default 0)
#' @param imax index of a hill from which summation stops (default the rest of hills)
#'
#' @export
#' @examples
#' tfes<-fes(acealanme)
#' minima<-fesminima(tfes)
#' prof<-feprof(minima)
#' prof
feprof <- function(minims=minims, imin=0, imax=NULL) {
  fes<-minims$fes
  rows<-minims$rows
  mins<-minims$minima
  hills<-minims$hills
  if(is.null(imax)) {
    imax<-nrow(hills)
    cat(imax,"\n")
  }
  if(imax>nrow(hills)) {
    imax<-nrow(hills)
    cat("Warning: You requested more hills by imax than available, using all hills\n")
  }
  if(imin>=imax) {
    stop("Error: imax must be higher than imin")
  }
  tt <- imin:imax
  mms <- data.frame(tt)
  if(minims$dimension==1) {
    for(i in 1:nrow(mins)) {
      if(minims$per[1]==T) {
        mm<-fe1dp(hills[,2], hills[,3], hills[,4], mins[i,3], minims$pcv1[2]-minims$pcv1[1], imin, imax)
      } else {
        mm<-fe1d(hills[,2], hills[,3], hills[,4], mins[i,3], imin, imax)
      }
      mms<-cbind(mms,mm)
    }
  }
  if(minims$dimension==2) {
    for(i in 1:nrow(mins)) {
      if(minims$per[1]==T && minims$per[2]==T) {
        mm<-fe2dp12(hills[,2], hills[,3], hills[,4], hills[,5], hills[,6], mins[i,4], mins[i,5], minims$pcv1[2]-minims$pcv1[1], minims$pcv2[2]-minims$pcv2[1], imin, imax)
      }
      if(minims$per[1]==T && minims$per[2]==F) {
        mm<-fe2dp1(hills[,2], hills[,3], hills[,4], hills[,5], hills[,6], mins[i,4], mins[i,5], minims$pcv1[2]-minims$pcv1[1], imin, imax)
      }
      if(minims$per[1]==F && minims$per[2]==T) {
        mm<-fe2dp2(hills[,2], hills[,3], hills[,4], hills[,5], hills[,6], mins[i,4], mins[i,5], minims$pcv2[2]-minims$pcv2[1], imin, imax)
      }
      if(minims$per[1]==F && minims$per[2]==F) {
        mm<-fe2d(hills[,2], hills[,3], hills[,4], hills[,5], hills[,6], mins[i,4], mins[i,5], imin, imax)
      }
      mms<-cbind(mms,mm)
    }
  }
  profs<-list(mms=mms, mins=mins, fes=fes, rows=rows, dimension=minims$dimension, per=minims$per, pcv1=minims$pcv1, pcv2=minims$pcv2)
  class(profs) <- "profiles"
  return(profs)
}

#' @export
print.profiles <- function(profs=profs) {
  cat(nrow(profs$mms))
  if(profs$dimension==1) {
    cat(" 1D minima\n")
  } else {
    cat(" 2D minima\n")
  }
}

#' Print summary for free energy profile
#'
#' `summary.profiles` prints the list of free energy minima with maximal
#' and minimal free energy differences.
#'
#' @param profs profiles object
#' @param imind index of a hill from which calculation of difference starts (default 0)
#' @param imaxd index of a hill from which calculation of difference stops (default the rest of hills)
#'
#' @export
#' @examples
#' tfes<-fes(acealanme)
#' minima<-fesminima(tfes)
#' prof<-feprof(minima)
#' summary(prof)
summary.profiles <- function(profs=profs, imind=1, imaxd=NULL) {
  if(is.null(imaxd)) {
    imaxd<-nrow(profs$mms)
  }
  if(profs$dimension==1) {
    outprofile <- profs$mins
    mms<-profs$mms[,2:ncol(profs$mms)]-profs$mms[,2]
    outprofile <- cbind(outprofile,apply(mms[imind:imaxd,],2,min))
    outprofile <- cbind(outprofile,apply(mms[imind:imaxd,],2,max))
    outprofile <- cbind(outprofile,t(mms[imaxd,]))
    names(outprofile)[5:7]<-c("min diff", "max diff", "tail")
    print(outprofile)
  } else {
    outprofile <- profs$mins
    mms<-profs$mms[,2:ncol(profs$mms)]-profs$mms[,2]
    outprofile <- cbind(outprofile,apply(mms[imind:imaxd,],2,min))
    outprofile <- cbind(outprofile,apply(mms[imind:imaxd,],2,max))
    outprofile <- cbind(outprofile,t(mms[imaxd,]))
    names(outprofile)[7:9]<-c("min diff", "max diff", "tail")
    print(outprofile)
  }
}

#' Plot free energy profile
#'
#' `plot.profiles` plots evolution of free energy differences between minima.
#' They are colored by rainbow colors from the global one (blue) to the highest (red).
#'
#' @param profs profiles object
#' @param which vector of indexes of profiles to be ploted (default all)
#' @param ignoretime time in the first column of the HILLS file will be ignored
#' @inherit plot
#' @inherit image
#' @inherit contours
#'
#' @export
#' @examples
#' tfes<-fes(acealanme)
#' minima<-fesminima(tfes)
#' prof<-feprof(minima)
#' plot(prof)
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

