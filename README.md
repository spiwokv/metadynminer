[![Build Status](https://github.com/spiwokv/metadynminer/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/spiwokv/metadynminer/actions/) 
[![Build Status](https://ci.appveyor.com/api/projects/status/github/spiwokv/metadynminer?branch=master&svg=true)](https://ci.appveyor.com/project/spiwokv/metadynminer) 
[![CRAN status](https://www.r-pkg.org/badges/version/metadynminer)](https://cran.r-project.org/package=metadynminer) 
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/metadynminer)](https://cran.r-project.org/package=metadynminer)
[![Rdoc](https://api.rdocumentation.org/badges/version/metadynminer)](https://www.rdocumentation.org/packages/metadynminer)
[![codecov](https://codecov.io/gh/spiwokv/metadynminer/branch/master/graph/badge.svg)](https://app.codecov.io/gh/spiwokv/metadynminer/)

# MetadynMiner

## Web site
https://metadynamics.cz/metadynminer/

https://spiwokv.github.io/metadynminer/ ([pkgdown](https://pkgdown.r-lib.org/) web)

## Introduction
MetadynMiner is R packages for reading, analysis and visualization of metadynamics HILLS files produced by Plumed.
It reads HILLS files from Plumed, calculates free energy surface by fast Bias Sum algorithm, finds minima and analyses
transition paths by Nudged Elastic Band method.

## Usage
```R
# Install from R repository
install.packages("metadynminer")

# Install from GitHub by devtools
install.packages("devtools")
devtools::install_github("spiwokv/metadynminer")

# Load library
library(metadynminer)
# Read hills file
hillsf<-read.hills("HILLS", per=c(TRUE, TRUE)) # HILLS with periodicity on CV1 and CV2

# Sum two hills files
hillsf+hillsf

# Summary of a hills file
summary(hillsf)

# Plot CVs
plot(hillsf)
```
![hills2d](./figs/hills2d.png)
```R
# Plot heights
plotheights(hillsf)
```
![hills2dh](./figs/hills2dh.png)
```R
# Calculate FES by bias sum (alternatively use fes2 for conventional calculation)
tfes<-fes(hillsf)

# Calculate FES for given range (indexes of hills)
tfes<-fes(hillsf, imin=5000, imax=10000)

# Sum two FESes
tfes+tfes

# Calculate and subtract min, max or mean from a FES
tfes<-tfes-min(tfes)

# Summary of FES
summary(tfes)

# Plot FES
plot(tfes)
```
![fes2d](./figs/fes2d.png)
```R
# Plot FES with color scale
plot(tfes, colscale=T)
```
![fes2d3](./figs/fes2d3.png)
```R
# Find minima
minima<-fesminima(tfes)

# Summary of minima
summary(minima)

# Create empty minima list, create ad hoc minimum, add minima
minima<-oneminimum(tfes, cv1=0, cv2=0)
fesminima(fes) + oneminimum(tfes, cv1=0, cv2=0)

# Plot free energy minima
plot(minima)
```
![min2d](./figs/min2d.png)
```R
# Calculate free energy profile for minima
prof<-feprof(minima)

# Plot free energy profile for minima
plot(prof)
```
![prof2d](./figs/prof2d.png)
```R

# Make 1D free energy surface from the 2D one
tfes1<-fes2d21d(hillsf, remdim=2) # T=300K, kJ/mol
plot(tfes1)
```
![fes1d](./figs/fes1d.png)
```R
# Calculate transition path using Nudged Elastic Band
myneb <- neb(minima, min1="A", min2="B")
myneb
summary(myneb)

# Plot transition path
plot(myneb)
```
![neb1](./figs/neb1.png)
```R
# Plot transition path on FES
plot(minima)
linesonfes(myneb)
```
![neb2](./figs/neb2.png)

## Tips and Tricks
### Publication quality figures
Following script can be used to generate a publication quality figure (8x8 cm, 600 dpi):
```R
hillsf <- read.hills("HILLS", per=c(T,T))
tfes<-fes(hillsf)
png("filename.png", height=8, width=8, units='cm', res=600, pointsize=6)
plot(tfes)
dev.off()
```
[Image](./figs/fes2dq.png)

### Making FES relative to the global minimum
You can set free energy minimum to zero by typing:
```R
tfes <- tfes - min(tfes)
```

### Hills from restarted simulations
There are two ways to cope with hills file from restarted simulations, i.e. simulations
where the time column starts from zero at every restart. It is possible to set `ignoretime=T`
in the `read.hills` function. It will take the time of the first hill and use it as the uniform
step.

Alternatively, it is possible to use `ignoretime=T` for other ploting functions.

MetadynMiner does not support hills files with variable PACE.

### Making movie
Individual snapshots of a movie can be generated and concatenated by:
```R
#install.packages("magick") <- install before the first run
library(metadynminer)
library(magick)
odir<-file.path(tempdir(), "frames")
dir.create(odir, recursive=T)
hillsf <- read.hills("HILLS", per=c(T,T))
tfes<-fes(hillsf, imax=300)
png(paste(odir, "/snap%04d.png", sep=""))
plot(tfes, zlim=c(-200,0))
for(i in 1:99) {
  tfes<-tfes+fes(hillsf, imin=300*i+1, imax=300*(i+1))
  plot(tfes, zlim=c(-200,0))
}
dev.off()
figs <- list.files(odir, full.names=T)
rfigs <- lapply(figs, image_read)
allfigs <- image_join(rfigs)
anim <- image_animate(allfigs, fps=25)
image_write(image=anim, path="anim.gif")
```
![anim1](./figs/anim.gif)

Package [magick](https://github.com/ropensci/magick) was used. Alternatively,
files can be concatenated by outside R by any movie making program.

If you instead want to see flooding, type:
```R
#install.packages("magick") <- install before the first run
library(metadynminer)
library(magick)
odir<-file.path(tempdir(), "frames")
dir.create(odir, recursive=T)
hillsf <- read.hills("HILLS", per=c(T,T))
tfes<-fes(hillsf)
png(paste(odir, "/snap%04d.png", sep=""))
plot(tfes, zlim=c(-200,0))
for(i in 0:99) {
  tfes<-tfes + -1*fes(hillsf, imin=300*i+1, imax=300*(i+1))
  plot(tfes, zlim=c(-200,0))
}
dev.off()
figs <- list.files(odir, full.names=T)
rfigs <- lapply(figs, image_read)
allfigs <- image_join(rfigs)
anim <- image_animate(allfigs, fps=25)
image_write(image=anim, path="flooding.gif")
```
![anim2](./figs/flooding.gif)


### Evaluation of convergence of one CV
You can use function `fes2d21d` to convert a 2D surface to 1D and to evaluate the evolution:
```R
hillsf <- read.hills("HILLS", per=c(T,T))
tfes1<-fes2d21d(hillsf, remdim=2)
plot(tfes1-min(tfes1), ylim=c(0,80), lwd=4, col="black")
for(i in 1:10) {
 tfes1<-fes2d21d(hillsf, imax=3000*i)
 lines(tfes1-min(tfes1), col=rainbow(13)[i])
}
```
![evol1d](./figs/evol1d.png)

### Transforming CVs
If you want to use degrees instead of radians on axes, set `axes=F` in the plot function and then plot
(without closing the plot window!) both axes separately.
```R
plot(tfes, axes=F)
axis(2, at=-3:3*pi/3, labels=-3:3*60)
axis(1, at=-3:3*pi/3, labels=-3:3*60)
box()
```
![degs](./figs/degs.png)

The expression -3:3 will generate a vector {-3,-2,-1,0,1,2,3}, which can be multiplied by pi/3
(tick positions in radians) or by 60 (tick positions in degrees). If you want to transform just
one axis, e.g. the horizontal one while keeping the vertical unchanged, simply type `axis(2)`
for the vertical one. `box()` redraws a box.

### kcal vs kJ
MetadynMiner works in kJ/mol by defauls. If your MD engine uses kcal/mol instead, you can either
multiply your free energy surface by 4.184 to get kJ/mol. If you prefer to keep kcal/mol, you can
set `eunit="kcal/mol"` for functions `fes2d21d` or `summary` of minima object. Other units are not
supported.

### Shifting a periodic CV
It may happen that some simulations with a torsion CV it may be difficult to analyze and visualize
it in the range -pi - +pi. However, this problem is not very common so we did not make any user
friendly way how to solve this and it can be solved in a user unfriendly way. Let us consider
we want to shift the first collective variable to be in the range 0 - 2pi. First we will make 
copy of acealanme. We will change its `pcv1` to `c(0,2*pi)`. Finally we can add 2pi to the
first collective variable:
```R
acealanmec<-acealanme
acealanmec$pcv1<-c(0,2*pi)
acealanmec$hillsfile[acealanmec$hillsfile[,2]<0,2]<-
    acealanmec$hillsfile[acealanmec$hillsfile[,2]<0,2]+2*pi
tfes<-fes(acealanmec)
plot(tfes)
```
![shift](./figs/shift.png)

The hills file object has several instances including `hillsfile`, which contains the HILLS file,
and `pcv1` with collective variable periodicity. They can be printed by `$` operator. The expression
`acealanmec$hillsfile[,2]` prints all values of the first collective variable. The expression
`acealanmec$hillsfile[,2]<0` prints the same number of `TRUE` or `FALSE` values depending whether
the first collective variable is positive or negative. The expression
`acealanmec$hillsfile[acealanmec$hillsfile[,2]<0,2]` prints only negative values of the first
collective variable. They can be replaced by the same value + 2pi.

## Training
[MetadynMiner Webinar video](https://youtu.be/W8N-G8d0or4)

[Plumed Masterclass tutorial](https://www.plumed.org/doc-master/user-doc/html/masterclass-22-02.html)

[Plumed Masterclass video 1](https://youtu.be/T8a-kP6V3_g)

[Plumed Masterclass video 2](https://youtu.be/q1D39A_LQag)

## Contact
Vojtech Spiwok - spiwokv{youknowwhat}vscht.cz

To contribute, se [CONTRIBUTING.md](https://github.com/spiwokv/metadynminer/blob/master/CONTRIBUTING.md)
  

