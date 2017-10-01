source("metadynminer.R")
hillsf <- read.hills("HILLS_1.3_full", per=c(T,T))
#hillsf <- read.hills("HILLS_AceProProNH2", per=c(F))
tfes<-fes2d(hillsf, tmax=1000)
png("test%03d.png")
plot(tfes, zlim=c(-120,0))
for(i in 1:10) {
  tfes<-tfes+fes2d(hillsf, tmin=1000*i, tmax=1000*(i+1))
  plot(tfes, zlim=c(-120,0))
}
#image(x=180*tfes$x/pi, y=180*tfes$y/pi, z=tfes$fes, axes=F)
#contour(x=180*tfes$x/pi, y=180*tfes$y/pi, z=tfes$fes, add=T)
#axis(1, at=-3:3*60)
#axis(2, at=-3:3*60)
#box()
dev.off()


