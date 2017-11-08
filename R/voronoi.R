x<-rnorm(10)
y<-rnorm(10)
triangles<-c(-5,-4,5,-4,0,5)
png("voronoi.png")
plot(x,y, xlim=c(-5,5), ylim=c(-5,5))
lines(triangles[c(0:2*2+1,1)], triangles[c(0:2*2+2,2)])
for(i in 1:10) {
  for(j in 1:nrow(triangles)) {
    
  }
}
dev.off()

