# MetadynMiner

## Introduction
MetadynMiner is R packages for reading, analysis and visualisation of metadynamics HILLS files produced by Plumed.

## Usage
```R
# Install from R repository
install.packages("metadynminer") # in future will be added to R repository

# Load library
library(metadynminer) # in future will be added to R repository

# Read hills file
hillsf<-read.hills("HILLS", per=c(TRUE, TRUE)) # HILLS with periodicity on CV1 and CV2

# Sum two hills files
hillsf+hillsf

# Summary of a hills file
summary(hillsf)

# Plot CVs
plot(hillsf)

# Plot heights
plotheights(hillsf)

# Calculate FES by bias sum (alternatively use fes2 for conventional calculation)
tfes<-fes(hillsf)

# Calculate FES for given range (indexes of hills)
tfes<-fes(hillsf, imin=5000, imax=10000)

# Sum two FESes
tfes+tfes

# Calculate and substract min, max or mean from a FES
tfes<-tfes-min(tfes)

# Summary of FES
summary(tfes)

# Plot FES
plot(tfes)

# Plot FES with color scale
plot(tfes, colscale=T)

# Find minima
minima<-fesminima(tfes)

# Summary of minima
summary(minima)

# Create empty minima list, create ad hoc minimum, add minima
minima<-emptyminima(tfes)
minima<-oneminimum(tfes, cv1=0, cv2=0)
fesminima(fes) + oneminimum(tfes, cv1=0, cv2=0)

# Plot free energy minima
plot(minima)

# Calculate free energy profile for minima
prof<-feprof(minima)

# Plot free energy profile for minima
plot(prof)

# Make 1D free energy surface from the 2D one
tfes1<-fes2d21d(hillsf, remdim=2) # T=300K, kJ/mol
plot(tfes1)

# Calculate transition path using Nudged Elastic Band
myneb <- neb(minima, min1="A", min2="B")
myneb
summary(myneb)

# Plot transition path
plot(myneb)

# Plot transition path on FES
plot(minima)
linesonfes(myneb)
```

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
Individual snapshots of a movie can be generated by:
```R
hillsf <- read.hills("HILLS", per=c(T,T))
tfes<-fes(hillsf, tmax=100)
png("snap%04d.png")
plot(tfes, zlim=c(-200,0))
for(i in 1:299) {
 tfes<-tfes+fes(hillsf, imin=100*i+1, imax=100*(i+1))
 plot(tfes, zlim=c(-200,0))
}
dev.off()
```
These files can be concatenated by a movie making program such as mencoder.

If you instead want to see flooding, type:
```R
hillsf <- read.hills("HILLS", per=c(T,T))
tfes<-fes(hillsf)
png("snap%04d.png")
plot(tfes, zlim=c(-200,0))
for(i in 0:299) {
  tfes<-tfes + -1*fes(hillsf, imin=100*i+1, imax=100*(i+1))
  plot(tfes, zlim=c(-200,0))
}
dev.off()
```

### Evaluation of convergence of one CV
You can use function `fes2d21d` to convert a 2D surface to 1D and to evalulate the evolution:
```R
hillsf <- read.hills("HILLS", per=c(T,T))
tfes1<-fes2d21d(hillsf, remdim=2)
plot(tfes1-min(tfes1), ylim=c(0,80), lwd=4, col="black")
for(i in 1:10) {
 tfes1<-fes2d21d(hillsf, imax=3000*i)
 lines(tfes1-min(tfes1), col=rainbow(13)[i])
}
```

### Transforming CVs
If you want to use degrees instead of radians on axes, set axes=F in the plot function and then plot
(without closing the plot window!) both axes separately.
```R
plot(tfes, axes=F)
axis(2, at=-3:3*pi/3, labels=-3:3*60)
axis(1, at=-3:3*pi/3, labels=-3:3*60)
```
The expression -3:3 will generate a vector {-3,-2,-1,0,1,2,3}, which can be multiplied by pi/3
(tick positions in radians) or by 60 (tick positions in degrees). If you want to transform just
one axis, e.g. the horizontal one while keepinng the vertical unchanges, simply type `axis(2)`
for the vertical one.

### kcal vs kJ
MetadynMiner works in kJ/mol by defauls. If your MD engine uses kcal/mol instead, you can either
multiply your free energy surface by 4.184 to get kJ/mol. If you prefer to keep kcal/mol, you can
set `eunit="kcal/mol"` for functions `fes2d21d` or `summary` of minima object. Other units are not
supported.

### Shifting a periodic CV
It may happen that some simulations with a torsion CV it may be difficult to analyse and visualize
it in the range -pi - +pi. However, this problem is not very common so we did not make any user friendly way
how to solve this and it can be solved in a user unfriendly way. It is possible to move the whole free
energy surface and corresponding x or y values within the free energy surface object. For example,
if you want to shift free energy surface to have phi from 0 to 2pi you can do this:
```R
hillsf <- read.hills("HILLS", per=c(T,T))
tfes<-fes(hillsf)
tfes$fes<- rbind(tfes$fes[129:256,], tfes$fes[1:128,])
  # replace 128, 129 and 256 if tfes$rows!=256
tfes2$x<-tfes2$x+pi
plot(tfes2)
```
