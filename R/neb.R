
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

# plot FES
plot.fes<-function(inputfes, plottype="both",
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
    x<-inputfes$x
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
    x<-inputfes$x
    y<-inputfes$y
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

print.minima<-function(minims) {
  cat("$minima\n\n")
  print(minims$minima)
}

neb<-function(minims, min1, min2, nbins=20,
              nsteps=100, step=1.0, k=0.2) {
  fes<-minims$fes
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
  path<-data.frame(pathx,pathy)
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
       newpath[j+1,1] = newpath[j+1,1] + step*(f1p[1]+f2a[1]+fphi*f2p[1])
       newpath[j+1,2] = newpath[j+1,2] + step*(f1p[2]+f2a[2]+fphi*f2p[2])
    }
    path <- newpath
  }
  return(path)
}

fes<-read.fes("bias.txt")
minima<-fesminima(fes)
png("test.png")
plot(fes)
text(minima$minima[,4], minima$minima[,5], labels=minima$minima[,1])
path<-neb(minima, min1=minima$minima[1,], min2=minima$minima[2,])
lines(pi*(path-65)/64)
points(pi*(path-65)/64)
path<-neb(minima, min1=minima$minima[1,], min2=minima$minima[3,])
lines(pi*(path-65)/64)
points(pi*(path-65)/64)
path<-neb(minima, min1=minima$minima[2,], min2=minima$minima[5,])
lines(pi*(path-65)/64)
points(pi*(path-65)/64)
path<-neb(minima, min1=minima$minima[3,], min2=minima$minima[5,])
lines(pi*(path-65)/64)
points(pi*(path-65)/64)
path<-neb(minima, min1=minima$minima[3,], min2=minima$minima[5,])
lines(pi*(path-65)/64)
points(pi*(path-65)/64)
path<-neb(minima, min1=minima$minima[2,], min2=minima$minima[4,])
lines(pi*(path-65)/64)
points(pi*(path-65)/64)
dev.off()

