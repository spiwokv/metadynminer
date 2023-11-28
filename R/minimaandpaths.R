#' Find free energy minima in the fes object (generic function for 'metadynminer'
#' and 'metadynminer3d')
#'
#' `fesminima` finds free energy minima on 1D or 2D free energy surface.
#' The surface is divided by a 1D or 2D grid and minima are found for each
#' bin. Next the program determines whether the minimum of a bin is a local
#' minimum of the whole free energy surface. Free energy minima are labeled
#' constitutively by capital letters.
#'
#' @param inputfes fes object.
#' @param nbins number of bins for each CV (default 8).
#' @return minima object.
#'
#' @export fesminima
fesminima<-function(inputfes, nbins) {
  UseMethod("fesminima")
}

#' Find free energy minima in the fes object
#'
#' `fesminima.fes` finds free energy minima on 1D or 2D free energy surface.
#' The surface is divided by a 1D or 2D grid and minima are found for each
#' bin. Next the program determines whether the minimum of a bin is a local
#' minimum of the whole free energy surface. Free energy minima are labeled
#' constitutively by capital letters.
#'
#' @param inputfes fes object.
#' @param nbins number of bins for each CV (default 8).
#' @return minima object.
#'
#' @export
#' @examples
#' tfes<-fes(acealanme, imax=5000)
#' minima<-fesminima(tfes)
#' minima
fesminima.fes<-function(inputfes, nbins=8) {
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

#' Creates one ad hoc free energy minimum for a fes object (generic function
#' for 'metadynminer' and 'metadynminer3d')
#'
#' `oneminimum` creates an ad hoc free energy minimum on free energy surface.
#' This can be used to calculate free energy surface evolution at arbitrary
#' point of free energy surface.
#'
#' @param inputfes fes object.
#' @param cv1 the value of collective variable 1.
#' @param cv2 the value of collective variable 2.
#' @param cv3 the value of collective variable 3.
#' @return minima object.
#'
#' @export oneminimum
oneminimum<-function(inputfes, cv1, cv2, cv3) {
  UseMethod("oneminimum")
}

#' Creates one ad hoc free energy minimum for a fes object
#'
#' `oneminimum.fes` creates an ad hoc free energy minimum on free energy surface.
#' This can be used to calculate free energy surface evolution at arbitrary
#' point of free energy surface.
#'
#' @param inputfes fes object.
#' @param cv1 the value of collective variable 1.
#' @param cv2 the value of collective variable 2.
#' @param cv3 the value of collective variable 3.
#' @return minima object.
#'
#' @export
#' @examples
#' tfes<-fes(acealanme1d)
#' minima<-fesminima(tfes)
#' minima<-minima+oneminimum(tfes, cv1=0, cv2=0)
#' minima
oneminimum.fes<-function(inputfes, cv1, cv2, cv3) {
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
    minima<-list(minima=minima, hills=min1$hills, fes=min1$fes, rows=min1$rows, dimension=min1$dimension, per=min1$per, x=min1$x, pcv1=min1$pcv1, pcv2=min1$pcv2)
    class(minima) <- "minima"
  }
  return(minima)
}

#' Print minima object
#'
#' `print.minima` prints free energy minima (identifier, values of bins and collective variables and free energy).
#'
#' @param x minima object.
#' @param ... further arguments passed to or from other methods.
#'
#' @export
#' @examples
#' tfes<-fes(acealanme, imax=5000)
#' minima<-fesminima(tfes)
#' minima
print.minima<-function(x,...) {
  print(x$minima)
}

#' Print minima object summary
#'
#' `summary.minima` prints summary for free energy minima (identifier, values of bins and collective variables,
#' free energy and equilibrium populations).
#'
#' @param object minima object
#' @param temp temperature in Kelvins
#' @param eunit energy units (kJ/mol or kcal/mol, kJ/mol is default)
#' @param ... further arguments passed to or from other methods.
#'
#' @export
#' @examples
#' tfes<-fes(acealanme, imax=5000)
#' minima<-fesminima(tfes)
#' summary(minima)
summary.minima<-function(object, temp=300, eunit="kJ/mol",...) {
  minims<-object
  toprint <- minims$minima
  tind = 6
  if(minims$dimension==1) {
    tind = 4
  }
  if(eunit=="kJ/mol") {
    toprint<-cbind(toprint, exp(-1000*toprint[,tind]/8.314/temp))
  }
  if(eunit=="kcal/mol") {
    toprint<-cbind(toprint, exp(-1000*toprint[,tind]*4.184/8.314/temp))
  }
  sumpop<-sum(toprint[,tind+1])
  toprint<-cbind(toprint, 100*toprint[,tind+1]/sumpop)
  names(toprint)[tind+1]<-"relative_pop"
  names(toprint)[tind+2]<-"pop"
  print(toprint)
}

#' Plot minima object
#'
#' `plot.minima` plots free energy surface with minima. The free energy surface is plotted the same
#' way as by plot.fes with additional minima labels.
#'
#' @param x minima object.
#' @param plottype specifies whether 2D free energy surface will be plotted
#'        as image, contours or both (default "both").
#' @param colscale specifies whether color scale will be plotted (default False).
#' @param colscalelab color scale label (default "free energy").
#' @param main an overall title for the plot: see 'title'.
#' @param sub a sub title for the plot: see 'title'.
#' @param xlab a title for the x axis: see 'title'.
#' @param ylab a title for the y axis: see 'title'.
#' @param pch plotting 'character', i.e., symbol to use. See 'points'
#' @param bg background (fill) color for the open plot symbols given by
#'        'pch = 21:25'.
#' @param cex character (or symbol) expansion: a numerical vector. This
#'        works as a multiple of 'par("cex")'.
#' @param asp the y/x aspect ratio, see 'plot.window'.
#' @param col color of the free energy surface. For 1D surface it is the color
#'        of the line. For 2D it is a list of colors such as that generated by
#'        'rainbow', 'heat.colors', 'topo.colors', 'terrain.colors' or similar
#'        functions (default=rainbow(135)[100:1]).
#' @param xlim numeric vector of length 2, giving the x coordinates range.
#' @param ylim numeric vector of length 2, giving the y coordinates range.
#' @param zlim numeric vector of length 2, giving the z coordinates range.
#' @param nlevels number of contour levels desired if 'levels' is not
#'        supplied.
#' @param levels numeric vector of levels at which to draw contour lines.
#' @param labels a vector giving the labels for the contour lines.  If 'NULL'
#'        then the levels are used as labels, otherwise this is coerced
#'        by 'as.character'.
#' @param labcex 'cex' for contour labeling. This is an absolute size, not a
#'        multiple of 'par("cex")'.
#' @param drawlabels logical. Contours are labeled if 'TRUE'.
#' @param method character string specifying where the labels will be located.
#'        Possible values are '"simple"', '"edge"' and '"flattest"'
#'        (the default). See the 'Details' section.
#' @param lwd contour line width.
#' @param contcol contour color.
#' @param lty line type for the lines drawn.
#' @param axes a logical value indicating whether both axes should be drawn
#'        on the plot.
#' @param textcol color of minima labels.
#' @param ... further arguments passed to or from other methods.
#'
#' @export
#' @examples
#' tfes<-fes(acealanme, imax=5000)
#' minima<-fesminima(tfes)
#' plot(minima)
plot.minima <- function(x, plottype="both",
                  xlim=NULL, ylim=NULL, zlim=NULL,
                  colscale=F, colscalelab="free energy",
                  main=NULL, sub=NULL,
                  xlab=NULL, ylab=NULL,
                  nlevels=10, levels=NULL,
                  col=rainbow(135)[100:1],
                  labels=NULL, labcex=0.6, drawlabels=TRUE,
                  method="flattest", textcol="black",
                  pch=1, bg="red", cex=1,
                  contcol=par("fg"), lty=par("lty"), lwd=par("lwd"),
                  asp=NULL, axes=TRUE,...) {
  minims <- x
  fes<-minims$fes
  rows<-minims$rows
  minlabs<-minims$minima[,1]
  if(minims$dimension==1) {
    x<-minims$x
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
        main=main, sub=sub, asp=asp)
    text(minpoints, labels=minlabs, col=textcol, xlim=xlim, ylim=ylim, cex=cex)
  } else {
    minpoints<-minims$minima[,4:5]
    x<-minims$x
    y<-minims$y
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
    if(colscale) {
      layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths=c(4,1))
    }
    if(plottype=="points") {
      text(minpoints, labels=minlabs, col=textcol, xlim=xlim, ylim=ylim,
        xlab=xlab, ylab=ylab, axes=axes,
        pch=pch, bg=bg, cex=cex,
        main=main, sub=sub, asp=asp)
    }
    if(plottype=="image" || plottype=="both") {
      image(x,y,fes, zlim=zlim,
        col=col, xlim=xlim, ylim=ylim,
        xlab=xlab, ylab=ylab, axes=axes,
        main=main, sub=sub, asp=asp)
      text(minpoints, labels=minlabs, col=textcol, xlim=xlim, ylim=ylim, cex=cex)
    }
    if(plottype=="contour") {
      contour(x,y,fes, zlim=zlim,
              nlevels=nlevels, levels=levels,
              labels=labels, labcex=labcex, drawlabels=drawlabels,
              method=method, col=contcol, lty=lty, lwd=lwd,
              main=main, sub=sub, asp=asp)
      text(minpoints, labels=minlabs, col=textcol, xlim=xlim, ylim=ylim, cex=cex)
    }
    if(plottype=="both") {
      contour(x,y,fes, zlim=zlim,
              nlevels=nlevels, levels=levels,
              labels=labels, labcex=labcex, drawlabels=drawlabels,
              method=method, col=contcol, lty=lty, lwd=lwd, add=T)
    }
    if(colscale) {
      smat<-matrix(seq(from=zlim[1], to=zlim[2], length.out=100))
      image(c(0), seq(from=zlim[1], to=zlim[2], length.out=100),
            t(smat), zlim=zlim, col=col, xlab="", ylab=colscalelab, axes=F)
      axis(2, lty=lty, lwd=lwd)
      box(lwd=lwd)
      par(mfrow=c(1,1))
    }
  }
}

#' Calculate free energy profile for minima object (generic function for 'metadynminer'
#' and 'metadynminer3d')
#'
#' `feprof` calculates free energy profiles for free energy minima. It finds the global minimum
#' at the `imax` and calculates the evolution of free energies of a local vs. the global free energy
#' minimum. The free energy of the global minimum is constant (zero).
#'
#' @param minims minima object.
#' @param imax index of a hill from which summation stops (default the rest of hills).
#'
#' @export
feprof <- function(minims, imax) {
  UseMethod("feprof")
}

#' Calculate free energy profile for minima object
#'
#' `feprof.minima` calculates free energy profiles for free energy minima. It finds the global minimum
#' at the `imax` and calculates the evolution of free energies of a local vs. the global free energy
#' minimum. The free energy of the global minimum is constant (zero).
#'
#' @param minims minima object.
#' @param imax index of a hill from which summation stops (default the rest of hills).
#'
#' @export
#' @examples
#' tfes<-fes(acealanme, imax=5000)
#' minima<-fesminima(tfes)
#' prof<-feprof(minima)
#' prof
feprof.minima <- function(minims, imax=NULL) {
  fes<-minims$fes
  rows<-minims$rows
  mins<-minims$minima
  hills<-minims$hills
  if(is.null(imax)) {
    imax<-nrow(hills)
  }
  if(imax>nrow(hills)) {
    imax<-nrow(hills)
    cat("Warning: You requested more hills by imax than available, using all hills\n")
  }
  tt <- 1:imax
  mms <- data.frame(tt)
  if(minims$dimension==1) {
    for(i in 1:nrow(mins)) {
      if(minims$per[1]==T) {
        mm<-fe1dp(hills[,2], hills[,3], hills[,4], mins[i,3], minims$pcv1[2]-minims$pcv1[1], 0, imax-1)
      } else {
        mm<-fe1d(hills[,2], hills[,3], hills[,4], mins[i,3], 0, imax-1)
      }
      mms<-cbind(mms,mm)
    }
  }
  if(minims$dimension==2) {
    for(i in 1:nrow(mins)) {
      if(minims$per[1]==T && minims$per[2]==T) {
        mm<-fe2dp12(hills[,2], hills[,3], hills[,4], hills[,5], hills[,6], mins[i,4], mins[i,5], minims$pcv1[2]-minims$pcv1[1], minims$pcv2[2]-minims$pcv2[1], 0, imax-1)
      }
      if(minims$per[1]==T && minims$per[2]==F) {
        mm<-fe2dp1(hills[,2], hills[,3], hills[,4], hills[,5], hills[,6], mins[i,4], mins[i,5], minims$pcv1[2]-minims$pcv1[1], 0, imax-1)
      }
      if(minims$per[1]==F && minims$per[2]==T) {
        mm<-fe2dp2(hills[,2], hills[,3], hills[,4], hills[,5], hills[,6], mins[i,4], mins[i,5], minims$pcv2[2]-minims$pcv2[1], 0, imax-1)
      }
      if(minims$per[1]==F && minims$per[2]==F) {
        mm<-fe2d(hills[,2], hills[,3], hills[,4], hills[,5], hills[,6], mins[i,4], mins[i,5], 0, imax-1)
      }
      mms<-cbind(mms,mm)
    }
  }
  profs<-list(mms=mms, mins=mins, fes=fes, rows=rows, dimension=minims$dimension, per=minims$per, pcv1=minims$pcv1, pcv2=minims$pcv2)
  class(profs) <- "profiles"
  return(profs)
}

#' Print profiles object
#'
#' `print.profiles` prints free energy profile.
#'
#' @param x minima object.
#' @param ... further arguments passed to or from other methods.
#'
#' @export
#' @examples
#' tfes<-fes(acealanme, imax=5000)
#' minima<-fesminima(tfes)
#' prof<-feprof(minima)
#' prof
print.profiles <- function(x,...) {
  outprofile <- x$mins
  print(outprofile)
}

#' Print summary for free energy profile
#'
#' `summary.profiles` prints the list of free energy minima with maximal
#' and minimal free energy differences.
#'
#' @param object profiles object.
#' @param imind index of a hill from which calculation of difference
#'        starts (default 1).
#' @param imaxd index of a hill from which calculation of difference
#'        stops (default the rest of hills).
#' @param ... further arguments passed to or from other methods.
#'
#' @export
#' @examples
#' tfes<-fes(acealanme, imax=5000)
#' minima<-fesminima(tfes)
#' prof<-feprof(minima)
#' summary(prof)
summary.profiles <- function(object, imind=1, imaxd=NULL,...) {
  profs<-object
  if(!is.null(imaxd)) {
    if(nrow(profs$mms)<imaxd) {
      cat("Warning: You requested more hills by imaxd than available, using all hills\n")
      imaxd<-nrow(profs$mms)
    }
  }
  if(is.null(imaxd)) {
    imaxd<-nrow(profs$mms)
  }
  if(imind>imaxd) {
    stop("Error: imaxd cannot be lower than imind")
  }
  if(profs$dimension==1) {
    outprofile <- profs$mins
    mms<-profs$mms[,2:ncol(profs$mms)]-profs$mms[,2]
    outprofile <- cbind(outprofile,apply(mms[imind:imaxd,],2,min))
    outprofile <- cbind(outprofile,apply(mms[imind:imaxd,],2,max))
    outprofile <- cbind(outprofile,t(mms[imaxd,]))
    names(outprofile)[5:7]<-c("min diff", "max diff", "tail")
    print(outprofile)
  } else if(profs$dimension==2) {
    outprofile <- profs$mins
    mms<-profs$mms[,2:ncol(profs$mms)]-profs$mms[,2]
    outprofile <- cbind(outprofile,apply(mms[imind:imaxd,],2,min))
    outprofile <- cbind(outprofile,apply(mms[imind:imaxd,],2,max))
    outprofile <- cbind(outprofile,t(mms[imaxd,]))
    names(outprofile)[7:9]<-c("min diff", "max diff", "tail")
    print(outprofile)
  } else if(profs$dimension==3) {
    outprofile <- profs$mins
    mms<-profs$mms[,2:ncol(profs$mms)]-profs$mms[,2]
    outprofile <- cbind(outprofile,apply(mms[imind:imaxd,],2,min))
    outprofile <- cbind(outprofile,apply(mms[imind:imaxd,],2,max))
    outprofile <- cbind(outprofile,t(mms[imaxd,]))
    names(outprofile)[9:11]<-c("min diff", "max diff", "tail")
    print(outprofile)
  }
}

#' Plot free energy profile
#'
#' `plot.profiles` plots evolution of free energy differences between minima.
#' They are colored by rainbow colors from the global one (blue) to the highest (red).
#'
#' @param x profiles object.
#' @param which vector of indexes of profiles to be plotted (default all).
#' @param ignoretime time in the first column of the HILLS file will be ignored.
#' @param main an overall title for the plot: see 'title'.
#' @param sub a sub title for the plot: see 'title'.
#' @param xlab a title for the x axis: see 'title'.
#' @param ylab a title for the y axis: see 'title'.
#' @param asp the y/x aspect ratio, see 'plot.window'.
#' @param col color code or name, see 'par'.
#' @param xlim numeric vector of length 2, giving the x coordinates range.
#' @param ylim numeric vector of length 2, giving the y coordinates range.
#' @param lwd line width.
#' @param axes a logical value indicating whether both axes should be drawn
#'        on the plot.
#' @param ... further arguments passed to or from other methods.
#'
#' @export
#' @examples
#' tfes<-fes(acealanme, imax=5000)
#' minima<-fesminima(tfes)
#' prof<-feprof(minima)
#' plot(prof)
plot.profiles <- function(x, which=NULL,
                          ignoretime=FALSE,
                          xlim=NULL, ylim=NULL,
                          main=NULL, sub=NULL,
                          xlab=NULL, ylab=NULL,
                          col=NULL, asp=NULL, lwd=1, axes=T,...) {
  profs<-x
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

#' Calculate free energy at given point in the CV space (generic function for 'metadynminer'
#' and 'metadynminer3d')
#'
#' `fespoint` calculates free energy at given point in the CV space 'coord'.
#' Hills are summed from 'imin' to `imax`. Printed output can be suppressed by setting
#' 'verb' to TRUE.
#'
#' @param hills hillsfile object.
#' @param coord coordinates of the point in the CV space.
#' @param imin index of a hill from which calculation of difference
#'        starts (default 1).
#' @param imax index of a hill from which summation stops (default the rest of hills).
#' @param verb if TRUE, the output is verbose (default TRUE).
#'
#' @export
fespoint <- function(hills, coord, imin, imax, verb) {
  UseMethod("fespoint")
}

#' Calculate free energy at given point in the CV space
#'
#' `fespoint.hillsfile` calculates free energy at given point in the CV space 'coord'.
#' Hills are summed from 'imin' to `imax`. Printed output can be suppressed by setting
#' 'verb' to TRUE.
#'
#' @param hills hillsfile object.
#' @param coord coordinates of the point in the CV space.
#' @param imin index of a hill from which calculation of difference
#'        starts (default 1).
#' @param imax index of a hill from which summation stops (default the rest of hills).
#' @param verb if TRUE, the output is verbose (default TRUE).
#'
#' @export
#' @examples
#' fespoint(acealanme, c(0,0), imax=5000)
fespoint.hillsfile <- function(hills, coord=NULL, imin=1, imax=NULL, verb=T) {
  if(!is.null(imax)) {
    if(hills$size[1]<imax) {
      cat("Warning: You requested more hills by imax than available, using all hills\n")
      imax<-hills$size[1]
    }
  }
  if(is.null(imax)) {
    imax<-hills$size[1]
  }
  if(imin==0) {
    imin <- 1
  }
  if(imax>0 & imin>imax) {
    stop("Error: imax cannot be lower than imin")
  }
  if(hills$size[2]==7) {
    if(length(coord)!=2) {
      stop("Error: for 2D fes you must use 2D vector as coord")
    }
    pcv1 <- hills$pcv1[2] - hills$pcv1[1]
    pcv2 <- hills$pcv2[2] - hills$pcv2[1]
    if((hills$per[1]==F)&(hills$per[2]==F)) {
      fe<-f2d(hills$hills[,2], hills$hills[,3], hills$hills[,4], hills$hills[,5],
                 hills$hills[,6], coord[1], coord[2], imin-1, imax-1)
    }
    if((hills$per[1]==T)&(hills$per[2]==F)) {
      if((coord[1]-max(hills$hills[,2]))>0.3*pcv1) {
        cat("Warning: coord quite outside periodic CV range")
      }
      if((min(hills$hills[,2])-coord[1])>0.3*pcv1) {
        cat("Warning: coord quite outside periodic CV range")
      }
      fe<-f2dp1(hills$hills[,2], hills$hills[,3], hills$hills[,4], hills$hills[,5],
                 hills$hills[,6], coord[1], coord[2], pcv1,
                 imin-1, imax-1)
    }
    if((hills$per[1]==F)&(hills$per[2]==T)) {
      if((coord[2]-max(hills$hills[,3]))>0.3*pcv2) {
        cat("Warning: coord quite outside periodic CV range")
      }
      if((min(hills$hills[,3])-coord[2])>0.3*pcv2) {
        cat("Warning: coord quite outside periodic CV range")
      }
      fe<-f2dp2(hills$hills[,2], hills$hills[,3], hills$hills[,4], hills$hills[,5],
                 hills$hills[,6], coord[1], coord[2], pcv2,
                 imin-1, imax-1)
    }
    if((hills$per[1]==T)&(hills$per[2]==T)) {
      if((coord[1]-max(hills$hills[,2]))>0.3*pcv1) {
        cat("Warning: coord quite outside periodic CV range")
      }
      if((min(hills$hills[,2])-coord[1])>0.3*pcv1) {
        cat("Warning: coord quite outside periodic CV range")
      }
      if((coord[2]-max(hills$hills[,3]))>0.3*pcv2) {
        cat("Warning: coord quite outside periodic CV range")
      }
      if((min(hills$hills[,3])-coord[2])>0.3*pcv2) {
        cat("Warning: coord quite outside periodic CV range")
      }
      fe<-f2dp12(hills$hills[,2], hills$hills[,3], hills$hills[,4], hills$hills[,5],
                  hills$hills[,6], coord[1], coord[2], pcv1, pcv2, imin-1, imax-1)
    }
    if(verb) {
      cat("free energy at the point: ")
      cat(coord[1])
      cat(" ")
      cat(coord[2])
      cat(" is ")
      cat(fe)
      cat("\n")
    }
  }
  if(hills$size[2]==5) {
    if(length(coord)!=1) {
      stop("Error: for 1D fes you must use 1D vector as coord")
    }
    pcv1 <- hills$pcv1[2] - hills$pcv1[1]
    if(hills$per[1]==F) {
      fe<-f1d(hills$hills[,2], hills$hills[,3], hills$hills[,4],
                coord[1], imin-1, imax-1)
    }
    if(hills$per[1]==T) {
      if((coord[1]-max(hills$hills[,2]))>0.3*pcv1) {
        cat("Warning: coord quite outside periodic CV range")
      }
      if((min(hills$hills[,2])-coord[1])>0.3*pcv1) {
        cat("Warning: coord quite outside periodic CV range")
      }
      fe<-f1dp(hills$hills[,2], hills$hills[,3], hills$hills[,4],
                coord[1], pcv1, imin-1, imax-1)
    }
    if(verb) {
      cat("free energy at the point: ")
      cat(coord[1])
      cat(" ")
      cat(coord[2])
      cat(" is ")
      cat(fe)
      cat("\n")
    }
  }
  return(fe)
}

