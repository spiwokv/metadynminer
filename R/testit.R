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
#png("fes2.png")
#plotheights(hillsf, ignoretime=T)
#plot(tfes)
mins<-fesminima(tfes)
print(mins)
#summary(mins)
#plot(mins)
png("profs.png")
prof1<-feprof(mins)
summary(prof1)
plot(prof1)
dev.off()




