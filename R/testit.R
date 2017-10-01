source("metadynminer.R")
hillsf <- read.hills("HILLS_1.3_full", per=c(F,F))
#hillsf <- read.hills("HILLS_AceProProNH2", per=c(F))
tfes<-fes2d(hillsf)
png("test.png")
#plot(tfes)
mins<-fesminima(tfes)
print(mins)
summary(mins)
plot(mins)
dev.off()


