# Metadynminer

Metadynminer is R packages for reading, analysis and visualisation of metadynamics HILLS and COLVAR files produced by Plumed.

```R
# Install from R repository
install.packages("metadynminer") # in future will be added to R repository

# Load library
library(metadynminer) # in future will be added to R repository

# Read hills file
hills<-read.hills("HILLS", per=c(TRUE, TRUE)) # done

# Sum two hills files
hills+hills # done

# Summary of a hills file
summary(hills) # done

# Plot CVs
plot(hills) # done

# Plot heights
plotheights(hills) # done

# Calculate FES by bias sum (alternatively use fes2 for conventional calculation)
fes<-fes(hills, tmin=5000, tmax=10000) # done, to do manual

# Sum two FESes
fes+fes # done, to do manual

# Calculate and substract min, max or mean from a FES
fes<-fes+min(fes1) # done, to do manual

# Summary of FES
summary(fes) # done, to do manual

# Plot FES
plot(fes) # done

# In order to make a movie you can use summation of FESes
tfes<-fes(hills, tmax=1000)
png("test%03d.png")
plot(tfes, zlim=c(-120,0))
for(i in 1:10) {
  tfes<-tfes+fes(hills, tmin=1000*i+1, tmax=1000*(i+1))
  plot(tfes, zlim=c(-120,0))
}
dev.off()

# Find minima
minima<-fesminima(fes) # done, to do manual

# Summary of minima
summary(minima) # done, to do manual

# Create empty minima
minima<-emptyminima(fes) # done, to do manual

# Create ad hoc minimum
minima<-oneminimum(fes, cv1=0, cv2=0) # done, to do manual

# sum minima
fesminima(fes) + oneminimum(fes, cv1=0, cv2=0) # done, to do manual

# Plot free energy minima
plot(minima) # 2D done, 1D to do, to do manual

# Calculate free energy profile for minima
prof<-feprofiles(hills, minima) # done, to do manual

# Find transition path
# Summary of transition path
# Minima areas
# Reweight
