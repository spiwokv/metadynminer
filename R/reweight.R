#' Calculate free energy surface from metadynamics by reweighting algorithm
#' by Tiwary and Parrinello, J. Phys. Chem. B (2015) <doi:10.1021/jp504920s>
#'
#' `reweighttiwary` calculates free energy surface from hills and colvars by
#' Tiwary and Parrinello algorithm.
#'
#' @param cvfile colvarfile object.
#' @param hillsfile hillsfile object.
#' @param nfes number of bias potential sums used to estimate c(t) (default 100).
#' @param temp temperature in Kelvins
#' @param eunit energy units (kJ/mol or kcal/mol, kJ/mol is default)
#' @param gamma bias factor of well-tempered metadynamics (default NULL).
#' @param imax index of a collective variable record from which summation stops (default the rest of hills).
#' @param xlim numeric vector of length 2, giving the CV1 coordinates range.
#' @param ylim numeric vector of length 2, giving the CV2 coordinates range.
#' @param npoints resolution of the free energy surface in number of points (default 60).
#' @param maxfe free energy of point on the output free energy surface with no sampling (default 100).
#' @param usefes2 logical parameter, if TRUE 'fes2' function is used instead of 'fes' (default FLASE).
#' @return fes object.
#'
#' @export
#' @examples
#' tfes<-reweighttiwary(cvs=acealanmeCVs, hills=acealanme, imax=5000)
reweighttiwary<-function(cvs, hills, npoints=60, maxfe=100, nfes=100, temp=300,
                         gamma=10, eunits="kJ/mol", imax=NULL, xlim=NULL, ylim=NULL, usefes2=F) {
  if(class(cvs)!="colvarfile") {
    stop("Error: Wrong colvarfile formate")
  }
  if(class(hills)!="hillsfile") {
    stop("Error: Wrong hillsfile formate")
  }
  if(cvs$ncvs!=(hills$size[2]-3)/2) {
    stop("Error: Different number of collective variables in colvarfile and hillsfile")
  }
  if(eunits=="kJ/mol") {
    beta <- 1000/8.314/temp
  } else {
    if(eunits=="kcal/mol") {
      beta <- 4184/8.314/temp
    } else {
      stop("Error: invalid energy unit, use kJ/mol or kcal/mol")
    }
  }
  if(cvs$ncvs==2) {
    minCV1 <- min(cvs$cvs[,1])
    maxCV1 <- max(cvs$cvs[,1])
    minCV2 <- min(cvs$cvs[,2])
    maxCV2 <- max(cvs$cvs[,2])
    xlims<-c(minCV1-0.05*(maxCV1-minCV1), maxCV1+0.05*(maxCV1-minCV1))
    ylims<-c(minCV2-0.05*(maxCV2-minCV2), maxCV2+0.05*(maxCV2-minCV2))
    if(!is.null(xlim)) {xlims<-xlim}
    if((hills$per[1]==T)&is.null(xlim)) {xlims<-hills$pcv1}
    if(!is.null(ylim)) {ylims<-ylim}
    if((hills$per[2]==T)&is.null(ylim)) {ylims<-hills$pcv2}
  }
  if(cvs$ncvs==1) {
    minCV1 <- min(cvs$cvs)
    maxCV1 <- max(cvs$cvs)
    xlims<-c(minCV1-0.05*(maxCV1-minCV1), maxCV1+0.05*(maxCV1-minCV1))
    if(!is.null(xlim)) {xlims<-xlim}
    if((hills$per[1]==T)&is.null(xlim)) {xlims<-hills$pcv1}
  }
  ebtac <- c()
  if(is.null(imax)) {
    ihills <- hills$size[1]
    framestosum <- hills$size[1]/nfes
  } else {
    ihills <- imax*hills$size[1]/length(cvs$times)
    framestosum <- 1:ihills/nfes
  }
  if(usefes2) {
    tfes <- fes2(hills, imax=1, npoints=npoints, xlim=xlims, ylim=ylims) -
            fes2(hills, imax=1, npoints=npoints, xlim=xlims, ylim=ylims)
  } else {
    suppressWarnings(tfes <- fes(hills, imax=1, npoints=npoints, xlim=xlims, ylim=ylims) -
                             fes(hills, imax=1, npoints=npoints, xlim=xlims, ylim=ylims))
  }
  for(i in 1:nfes) {
    if(usefes2) {
      tfes <- tfes + fes2(hills, imin=(i-1)*framestosum+1, imax=i*framestosum,
                          npoints=npoints, xlim=xlims, ylim=ylims)
    } else {
      tfes <- tfes + fes(hills, imin=(i-1)*framestosum+1, imax=i*framestosum,
                         npoints=npoints, xlim=xlims, ylim=ylims)
    }
    s1 <- sum(exp(-beta*tfes$fes))
    s2 <- sum(exp(-beta*tfes$fes/gamma))
    ebtac<-c(ebtac,s1/s2)
  }
  if(!is.null(cvs$bias)) {
    if(is.null(imax)) {
      bp <- cvs$bias
    } else {
      bp <- cvs$bias[1:imax]
    }
  } else {
    stop("Error: Input colvarfile does not contain bias potential")
  }
  if(cvs$ncvs>1) {
    nlines<-nrow(cvs$cvs)
  } else {
    nlines<-length(cvs$cvs)
  }
  ebtacc <- rep(ebtac, each=nlines/nfes)
  if(length(ebtacc)>nlines) ebtacc<-ebtacc[1:nlines]
  if(length(ebtacc)<nlines) ebtacc[length(ebtacc):nlines]<-ebtacc[length(ebtacc)]
  if(cvs$ncvs==2) {
    x<-0:(npoints-1)*(xlims[2]-xlims[1])/(npoints-1)+xlims[1]
    y<-0:(npoints-1)*(ylims[2]-ylims[1])/(npoints-1)+ylims[1]
    if(is.null(imax)) {
      icv1 <- ceiling((cvs$cvs[,1]-xlims[1])*npoints/(xlims[2]-xlims[1]))
      icv2 <- ceiling((cvs$cvs[,2]-ylims[1])*npoints/(ylims[2]-ylims[1]))
    } else {
      icv1 <- ceiling((cvs$cvs[1:imax,1]-xlims[1])*npoints/(xlims[2]-xlims[1]))
      icv2 <- ceiling((cvs$cvs[1:imax,2]-ylims[1])*npoints/(ylims[2]-ylims[1]))
    }
    ofes <- matrix(rep(0, npoints*npoints), nrow=npoints)
    for(i in 1:npoints) {
      for(j in 1:npoints) {
        ofes[i,j]<-sum((icv1==i)*(icv2==j)*exp(beta*bp)/ebtacc)
      }
    }
    ofes <- -log(ofes)/beta
    ofes <- ofes - min(ofes)
    ofes[ofes==Inf]<-maxfe
    cfes<-list(fes=ofes, hills=hills$hillsfile, rows=npoints, dimension=2, per=hills$per, x=x, y=y, pcv1=hills$pcv1, pcv2=hills$pcv2)
    class(cfes) <- "fes"
  }
  if(cvs$ncvs==1) {
    minCV1 <- min(cvs$cvs)
    maxCV1 <- max(cvs$cvs)
    xlims<-c(minCV1-0.05*(maxCV1-minCV1), maxCV1+0.05*(maxCV1-minCV1))
    if(!is.null(xlim)) {xlims<-xlim}
    if((hills$per[1]==T)&is.null(xlim)) {xlims<-hills$pcv1}
    x<-0:(npoints-1)*(xlims[2]-xlims[1])/(npoints-1)+xlims[1]
    if(is.null(imax)) {
      icv1 <- ceiling((cvs$cvs-xlims[1])*npoints/(xlims[2]-xlims[1]))
    } else {
      icv1 <- ceiling((cvs$cvs[1:imax]-xlims[1])*npoints/(xlims[2]-xlims[1]))
    }
    ofes <- rep(0, npoints)
    for(i in 1:npoints) {
      ofes[i]<-sum((icv1==i)*exp(beta*bp)/ebtacc)
    }
    ofes <- -log(ofes)/beta
    ofes <- ofes - min(ofes)
    ofes[ofes==Inf]<-maxfe
    cfes<-list(fes=ofes, hills=hills$hillsfile, rows=npoints, dimension=1, per=hills$per, x=x, pcv1=hills$pcv1)
    class(cfes) <- "fes"
  }
  return(cfes)
}

#' Calculate free energy surface from biased simulations by reweighting algorithm
#' by Bonomi et al., J. Comput. Chem. (2009) <doi:10.1002/jcc.21305>
#'
#' `reweightbonomi` calculates free energy surface from colvars by
#' Bonomi et al. algorithm.
#'
#' @param cvfile colvarfile object.
#' @param temp temperature in Kelvins
#' @param eunit energy units (kJ/mol or kcal/mol, kJ/mol is default)
#' @param imin index of a collective variable record from which summation starts (default 1).
#' @param imax index of a collective variable record from which summation stops (default the rest of colvar values).
#' @param xlim numeric vector of length 2, giving the CV1 coordinates range.
#' @param ylim numeric vector of length 2, giving the CV2 coordinates range.
#' @param npoints resolution of the free energy surface in number of points (default 60).
#' @param maxfe free energy of point on the output free energy surface with no sampling (default 100).
#' @return fes object.
#'
#' @export
#' @examples
#' tfes<-reweightbonomi(cvs=acealanmeCVs, imax=5000)
reweightbonomi<-function(cvs, npoints=60, maxfe=100, temp=300,
                         eunits="kJ/mol", imin=1, imax=NULL, xlim=NULL, ylim=NULL) {
  if(class(cvs)!="colvarfile") {
    stop("Error: Wrong colvarfile formate")
  }
  if(eunits=="kJ/mol") {
    beta <- 1000/8.314/temp
  } else {
    if(eunits=="kcal/mol") {
      beta <- 4184/8.314/temp
    } else {
      stop("Error: invalid energy unit, use kJ/mol or kcal/mol")
    }
  }
  if(cvs$ncvs==2) {
    minCV1 <- min(cvs$cvs[,1])
    maxCV1 <- max(cvs$cvs[,1])
    minCV2 <- min(cvs$cvs[,2])
    maxCV2 <- max(cvs$cvs[,2])
    xlims<-c(minCV1-0.05*(maxCV1-minCV1), maxCV1+0.05*(maxCV1-minCV1))
    ylims<-c(minCV2-0.05*(maxCV2-minCV2), maxCV2+0.05*(maxCV2-minCV2))
    if(!is.null(xlim)) {xlims<-xlim}
    if(!is.null(ylim)) {ylims<-ylim}
  }
  if(cvs$ncvs==1) {
    minCV1 <- min(cvs$cvs)
    maxCV1 <- max(cvs$cvs)
    xlims<-c(minCV1-0.05*(maxCV1-minCV1), maxCV1+0.05*(maxCV1-minCV1))
    if(!is.null(xlim)) {xlims<-xlim}
  }
  if(is.null(imin)) {
    iminc <- 1
  } else {
    iminc <- imin
  }
  if(cvs$ncvs>1) {
    nlines<-nrow(cvs$cvs)
  } else {
    nlines<-length(cvs$cvs)
  }
  if(is.null(imax)) {
    imaxc <- nlines
  } else {
    imaxc <- imax
  }
  if(iminc>imaxc) stop("Error: imin > imax")
  if(imaxc>nlines) stop("Error: imax too high")
  if(!is.null(cvs$bias)) {
    bp <- cvs$bias[iminc:imaxc]
  } else {
    stop("Error: Input colvarfile does not contain bias potential")
  }
  if(cvs$ncvs==2) {
    x<-0:(npoints-1)*(xlims[2]-xlims[1])/(npoints-1)+xlims[1]
    y<-0:(npoints-1)*(ylims[2]-ylims[1])/(npoints-1)+ylims[1]
    icv1 <- ceiling((cvs$cvs[iminc:imaxc,1]-xlims[1])*npoints/(xlims[2]-xlims[1]))
    icv2 <- ceiling((cvs$cvs[iminc:imaxc,2]-ylims[1])*npoints/(ylims[2]-ylims[1]))
    ofes <- matrix(rep(0, npoints*npoints), nrow=npoints)
    for(i in 1:npoints) {
      for(j in 1:npoints) {
        ofes[i,j]<-sum((icv1==i)*(icv2==j)*exp(beta*bp))
      }
    }
    ofes <- -log(ofes)/beta
    ofes <- ofes - min(ofes)
    ofes[ofes==Inf]<-maxfe
    cfes<-list(fes=ofes, hills=NULL, rows=npoints, dimension=2, per=c(F,F), x=x, y=y, pcv1=c(-pi,pi), pcv2=c(-pi,pi))
    class(cfes) <- "fes"
  }
  if(cvs$ncvs==1) {
    minCV1 <- min(cvs$cvs)
    maxCV1 <- max(cvs$cvs)
    xlims<-c(minCV1-0.05*(maxCV1-minCV1), maxCV1+0.05*(maxCV1-minCV1))
    if(!is.null(xlim)) {xlims<-xlim}
    x<-0:(npoints-1)*(xlims[2]-xlims[1])/(npoints-1)+xlims[1]
    icv1 <- ceiling((cvs$cvs[iminc:imaxc]-xlims[1])*npoints/(xlims[2]-xlims[1]))
    ofes <- rep(0, npoints)
    for(i in 1:npoints) {
      ofes[i]<-sum((icv1==i)*exp(beta*bp))
    }
    ofes <- -log(ofes)/beta
    ofes <- ofes - min(ofes)
    ofes[ofes==Inf]<-maxfe
    cfes<-list(fes=ofes, hills=NULL, rows=npoints, dimension=1, per=c(F), x=x, pcv1=c(-pi,pi))
    class(cfes) <- "fes"
  }
  return(cfes)
}

#' Calculate free energy surface from biased simulations by reweighting algorithm
#' by Torrie and Valleau, J. Comput. Phys. (1977) <doi:10.1016/0021-9991(77)90121-8>
#'
#' `reweighttv` calculates free energy surface from colvars by Umbrella Sampling
#' reweighting algorithm. It should be used for a biased simulation with a static
#' bias potential. For metadynamics it is better to use 'reweighttiwary' or
#' 'reweightbonomi'.
#'
#' @param cvfile colvarfile object.
#' @param temp temperature in Kelvins
#' @param eunit energy units (kJ/mol or kcal/mol, kJ/mol is default)
#' @param imin index of a collective variable record from which summation starts (default 1).
#' @param imax index of a collective variable record from which summation stops (default the rest of colvar values).
#' @param xlim numeric vector of length 2, giving the CV1 coordinates range.
#' @param ylim numeric vector of length 2, giving the CV2 coordinates range.
#' @param npoints resolution of the free energy surface in number of points (default 60).
#' @param maxfe free energy of point on the output free energy surface with no sampling (default 100).
#' @return fes object.
#'
#' @export
#' @examples
#' tfes<-reweighttv(cvs=acealanmeCVs, imax=5000)
reweighttv<-function(cvs, npoints=60, maxfe=100, temp=300,
                         eunits="kJ/mol", imin=1, imax=NULL, xlim=NULL, ylim=NULL) {
  cat("Warning: Umbrella sampling reweighting should be used with a static bias potential!\n")
  cat("For metadynamics it is better to use 'reweighttiwary' or 'reweightbonomi'.\n")
  if(class(cvs)!="colvarfile") {
    stop("Error: Wrong colvarfile formate")
  }
  if(eunits=="kJ/mol") {
    beta <- 1000/8.314/temp
  } else {
    if(eunits=="kcal/mol") {
      beta <- 4184/8.314/temp
    } else {
      stop("Error: invalid energy unit, use kJ/mol or kcal/mol")
    }
  }
  if(cvs$ncvs==2) {
    minCV1 <- min(cvs$cvs[,1])
    maxCV1 <- max(cvs$cvs[,1])
    minCV2 <- min(cvs$cvs[,2])
    maxCV2 <- max(cvs$cvs[,2])
    xlims<-c(minCV1-0.05*(maxCV1-minCV1), maxCV1+0.05*(maxCV1-minCV1))
    ylims<-c(minCV2-0.05*(maxCV2-minCV2), maxCV2+0.05*(maxCV2-minCV2))
    if(!is.null(xlim)) {xlims<-xlim}
    if(!is.null(ylim)) {ylims<-ylim}
  }
  if(cvs$ncvs==1) {
    minCV1 <- min(cvs$cvs)
    maxCV1 <- max(cvs$cvs)
    xlims<-c(minCV1-0.05*(maxCV1-minCV1), maxCV1+0.05*(maxCV1-minCV1))
    if(!is.null(xlim)) {xlims<-xlim}
  }
  if(is.null(imin)) {
    iminc <- 1
  } else {
    iminc <- imin
  }
  if(cvs$ncvs>1) {
    nlines<-nrow(cvs$cvs)
  } else {
    nlines<-length(cvs$cvs)
  }
  if(is.null(imax)) {
    imaxc <- nlines
  } else {
    imaxc <- imax
  }
  if(iminc>imaxc) stop("Error: imin > imax")
  if(imaxc>nlines) stop("Error: imax too high")
  if(!is.null(cvs$bias)) {
    bp <- cvs$bias[iminc:imaxc]
  } else {
    stop("Error: Input colvarfile does not contain bias potential")
  }
  if(cvs$ncvs==2) {
    x<-0:(npoints-1)*(xlims[2]-xlims[1])/(npoints-1)+xlims[1]
    y<-0:(npoints-1)*(ylims[2]-ylims[1])/(npoints-1)+ylims[1]
    icv1 <- ceiling((cvs$cvs[iminc:imaxc,1]-xlims[1])*npoints/(xlims[2]-xlims[1]))
    icv2 <- ceiling((cvs$cvs[iminc:imaxc,2]-ylims[1])*npoints/(ylims[2]-ylims[1]))
    ofes <- matrix(rep(0, npoints*npoints), nrow=npoints)
    for(i in 1:npoints) {
      for(j in 1:npoints) {
        ofes[i,j]<-sum((icv1==i)*(icv2==j)*exp(beta*bp))
      }
    }
    ofes <- -log(ofes)/beta
    ofes <- ofes - min(ofes)
    ofes[ofes==Inf]<-maxfe
    cfes<-list(fes=ofes, hills=NULL, rows=npoints, dimension=2, per=c(F,F), x=x, y=y, pcv1=c(-pi,pi), pcv2=c(-pi,pi))
    class(cfes) <- "fes"
  }
  if(cvs$ncvs==1) {
    minCV1 <- min(cvs$cvs)
    maxCV1 <- max(cvs$cvs)
    xlims<-c(minCV1-0.05*(maxCV1-minCV1), maxCV1+0.05*(maxCV1-minCV1))
    if(!is.null(xlim)) {xlims<-xlim}
    x<-0:(npoints-1)*(xlims[2]-xlims[1])/(npoints-1)+xlims[1]
    icv1 <- ceiling((cvs$cvs[iminc:imaxc]-xlims[1])*npoints/(xlims[2]-xlims[1]))
    ofes <- rep(0, npoints)
    for(i in 1:npoints) {
      ofes[i]<-sum((icv1==i)*exp(beta*bp))
    }
    ofes <- -log(ofes)/beta
    ofes <- ofes - min(ofes)
    ofes[ofes==Inf]<-maxfe
    cfes<-list(fes=ofes, hills=NULL, rows=npoints, dimension=1, per=c(F), x=x, pcv1=c(-pi,pi))
    class(cfes) <- "fes"
  }
  return(cfes)
}

#' Calculate free energy surface from unbiased simulations as -kT log(P(s)).
#'
#' `weightgibbs` calculates free energy surface from colvars as -kT log(P(s)).
#' It should be used for unbiased simulations only. For biased simulations with
#' a static bias potential use 'reweightbonomi'. For metadynamics it is better
#' to use 'reweighttiwary'.
#'
#' @param cvfile colvarfile object.
#' @param temp temperature in Kelvins
#' @param eunit energy units (kJ/mol or kcal/mol, kJ/mol is default)
#' @param imin index of a collective variable record from which summation starts (default 1).
#' @param imax index of a collective variable record from which summation stops (default the rest of colvar values).
#' @param xlim numeric vector of length 2, giving the CV1 coordinates range.
#' @param ylim numeric vector of length 2, giving the CV2 coordinates range.
#' @param npoints resolution of the free energy surface in number of points (default 60).
#' @param maxfe free energy of point on the output free energy surface with no sampling (default 100).
#' @return fes object.
#'
#' @export
#' @examples
#' tfes<-weightgibbs(cvs=acealanmeCVs, imax=5000)
weightgibbs<-function(cvs, npoints=60, maxfe=100, temp=300,
                      eunits="kJ/mol", imin=1, imax=NULL, xlim=NULL, ylim=NULL) {
  if(class(cvs)!="colvarfile") {
    stop("Error: Wrong colvarfile formate")
  }
  if(eunits=="kJ/mol") {
    beta <- 1000/8.314/temp
  } else {
    if(eunits=="kcal/mol") {
      beta <- 4184/8.314/temp
    } else {
      stop("Error: invalid energy unit, use kJ/mol or kcal/mol")
    }
  }
  if(cvs$ncvs==2) {
    minCV1 <- min(cvs$cvs[,1])
    maxCV1 <- max(cvs$cvs[,1])
    minCV2 <- min(cvs$cvs[,2])
    maxCV2 <- max(cvs$cvs[,2])
    xlims<-c(minCV1-0.05*(maxCV1-minCV1), maxCV1+0.05*(maxCV1-minCV1))
    ylims<-c(minCV2-0.05*(maxCV2-minCV2), maxCV2+0.05*(maxCV2-minCV2))
    if(!is.null(xlim)) {xlims<-xlim}
    if(!is.null(ylim)) {ylims<-ylim}
  }
  if(cvs$ncvs==1) {
    minCV1 <- min(cvs$cvs)
    maxCV1 <- max(cvs$cvs)
    xlims<-c(minCV1-0.05*(maxCV1-minCV1), maxCV1+0.05*(maxCV1-minCV1))
    if(!is.null(xlim)) {xlims<-xlim}
  }
  if(is.null(imin)) {
    iminc <- 1
  } else {
    iminc <- imin
  }
  if(cvs$ncvs>1) {
    nlines<-nrow(cvs$cvs)
  } else {
    nlines<-length(cvs$cvs)
  }
  if(is.null(imax)) {
    imaxc <- nlines
  } else {
    imaxc <- imax
  }
  if(iminc>imaxc) stop("Error: imin > imax")
  if(imaxc>nlines) stop("Error: imax too high")
  if(!is.null(cvs$bias)) {
    cat("Warning: There is a bias potential column in colvar file.\n")
    cat("Consider using reweighting.\n")
  }
  if(cvs$ncvs==2) {
    x<-0:(npoints-1)*(xlims[2]-xlims[1])/(npoints-1)+xlims[1]
    y<-0:(npoints-1)*(ylims[2]-ylims[1])/(npoints-1)+ylims[1]
    icv1 <- ceiling((cvs$cvs[iminc:imaxc,1]-xlims[1])*npoints/(xlims[2]-xlims[1]))
    icv2 <- ceiling((cvs$cvs[iminc:imaxc,2]-ylims[1])*npoints/(ylims[2]-ylims[1]))
    ofes <- matrix(rep(0, npoints*npoints), nrow=npoints)
    for(i in 1:npoints) {
      for(j in 1:npoints) {
        ofes[i,j]<-sum((icv1==i)*(icv2==j))
      }
    }
    ofes <- -log(ofes)/beta
    ofes <- ofes - min(ofes)
    ofes[ofes==Inf]<-maxfe
    cfes<-list(fes=ofes, hills=NULL, rows=npoints, dimension=2, per=c(F,F), x=x, y=y, pcv1=c(-pi,pi), pcv2=c(-pi,pi))
    class(cfes) <- "fes"
  }
  if(cvs$ncvs==1) {
    minCV1 <- min(cvs$cvs)
    maxCV1 <- max(cvs$cvs)
    xlims<-c(minCV1-0.05*(maxCV1-minCV1), maxCV1+0.05*(maxCV1-minCV1))
    if(!is.null(xlim)) {xlims<-xlim}
    x<-0:(npoints-1)*(xlims[2]-xlims[1])/(npoints-1)+xlims[1]
    icv1 <- ceiling((cvs$cvs[iminc:imaxc]-xlims[1])*npoints/(xlims[2]-xlims[1]))
    ofes <- rep(0, npoints)
    for(i in 1:npoints) {
      ofes[i]<-sum(icv1==i)
    }
    ofes <- -log(ofes)/beta
    ofes <- ofes - min(ofes)
    ofes[ofes==Inf]<-maxfe
    cfes<-list(fes=ofes, hills=NULL, rows=npoints, dimension=1, per=c(F), x=x, pcv1=c(-pi,pi))
    class(cfes) <- "fes"
  }
  return(cfes)
}

