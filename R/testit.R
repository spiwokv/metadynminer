source("metadynminer.R")
hillsf <- read.hills("HILLS_1.3_full", per=c(T,T))
tfes<-fes2d(hillsf)
print("ok")
png("test.png")
image(tfes)
dev.off()


