source("metadynminer.R")
hillsf <- read.hills("HILLS", per=c(T,T))
#hillsf <- read.hills("HILLS_AceProProNH2", per=c(T))
tfes<-fes(hillsf)
# In order to make a movie you can use summation of FESes
#tfes<-fes(hillsf, tmax=100)
#png("frame%03d.png")
#plot(tfes, zlim=c(-180,0))
#for(i in 1:300) {
#  tfes<-tfes+fes(hillsf, tmin=100*i+1, tmax=100*(i+1))
#  plot(tfes, zlim=c(-180,0))
#}
#dev.off()
png("fes2.png")
#plotheights(hillsf, ignoretime=T)
#plot(tfes)
mins<-fesminima(tfes)
library(deldir)
mincoords<-mins$minima[,4:5]
mincoords<-rbind(mincoords,mins$minima[,4:5]+c(-2*pi,-2*pi))
mincoords<-rbind(mincoords,mins$minima[,4:5]+c(-2*pi,    0))
mincoords<-rbind(mincoords,mins$minima[,4:5]+c(-2*pi,+2*pi))
mincoords<-rbind(mincoords,mins$minima[,4:5]+c(    0,-2*pi))
mincoords<-rbind(mincoords,mins$minima[,4:5]+c(    0,+2*pi))
mincoords<-rbind(mincoords,mins$minima[,4:5]+c(+2*pi,-2*pi))
mincoords<-rbind(mincoords,mins$minima[,4:5]+c(+2*pi,    0))
mincoords<-rbind(mincoords,mins$minima[,4:5]+c(+2*pi,+2*pi))
vtess <- deldir(mincoords[,1],mincoords[,2])
#summary(mins)
plot(mins)
#plot(vtess, wlines="tess", wpoints="none", add=T)
border<-c()
for(i in 1:nrow(vtess$dirsgs)) {
  #if(sum(vtess$dirsgs[i,5]==(0:7*nrow(mins$minima)+1))) {
    #border<-rbind(border,c(vtess$dirsgs[i,c(1,3)],vtess$dirsgs[i,c(2,4)]))
    lines(vtess$dirsgs[i,c(1,3)],vtess$dirsgs[i,c(2,4)], col="red", lwd=2)
  #}
  #if(sum(vtess$dirsgs[i,6]==(0:7*nrow(mins$minima)+1))) {
#    border<-rbind(border,c(vtess$dirsgs[i,c(1,3)],vtess$dirsgs[i,c(2,4)]))
#  }
}
for(i in 1:nrow(vtess$delsgs)) {
  if(vtess$delsgs[i,5]<=nrow(mins$minima)) {
    lines(vtess$delsgs[i,c(1,3)],vtess$delsgs[i,c(2,4)], col="blue", lwd=2)
    print(vtess$delsgs[i,5:6])
  }
  if(vtess$delsgs[i,6]<=nrow(mins$minima)) {
    lines(vtess$delsgs[i,c(1,3)],vtess$delsgs[i,c(2,4)], col="blue", lwd=2)
    print(vtess$delsgs[i,5:6])
  }
}
#png("profs.png")
#prof1<-feprof(mins)
#summary(prof1)
#plot(prof1)
#dev.off()




