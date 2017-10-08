source("metadynminer.R")
hillsf <- read.hills("HILLS_1.3_full", per=c(T,T))
#hillsf <- read.hills("HILLS_AceProProNH2", per=c(T))
tfes<-fes(hillsf)
png("test.png")
#plot(tfes)
mins<-fesminima(tfes)
#print(mins)
#summary(mins)
#plot(mins)
prof1<-feprof(mins)
plot(prof1[,4]-prof1[,2])
dev.off()




