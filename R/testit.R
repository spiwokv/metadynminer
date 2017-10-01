source("metadynminer.R")
hillsf <- read.hills("HILLS_1.3_full", per=c(T,T))
tfes<-fes2d(hillsf)
print("ok")
png("test.png")
image(x=180*tfes$x/pi, y=180*tfes$y/pi, z=tfes$fes, axes=F)
axis(1, at=-3:3*60)
axis(2, at=-3:3*60)
box()
dev.off()


