# Metadynminer

Metadynminer is R packages for reading, analysis and visualisation of metadynamics HILLS files produced by Plumed.

> install.packages("metadynminer")
> library(metadynminer)

<> # Reading hills file
> hills<-read.hills("HILLS", per=c(TRUE, TRUE)) # done

# Summary
> summary(hills) # done

# Plot hills file
> plot(hills) # done

# Plot hill heights
> plotheights(hills) # done

# Calculate FES
> fes<-mtdfes(hills, time=10000)

# Calculate FESes during flooding
> fes2<-mtdfes2(hills, time=1:10*1000)

# Evaluate FES
> summary(fes) # done

# Plot FES
> plot(fes) # done

# Find minima
> minima<-fesminima(fes) # done

# Evaluate free energy minima
> summary(minima) # done

# Plot free energy minima
> plot(fes) # done

# Calculate free energy
> fe<-freeene(x, hills, per=c(TRUE, TRUE), time=10000)

# Calculate free energy profile
> fe<-freeene2(x, hills, per=c(TRUE, TRUE), time=0:10000)

# Calculate free energy difference profile
> fe1<-freeene2(x1, hills, per=c(TRUE, TRUE), time=0:10000)
> fe2<-freeene2(x2, hills, per=c(TRUE, TRUE), time=0:10000)
> fediff<-fe2-fe1

# Check if converged
> isconverged(fediff, from=5000)

# Objects:
# hillsfile: HILLS file
# fes: free energy surface
# feminima: free energy minima

