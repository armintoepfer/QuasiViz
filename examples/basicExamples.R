setwd("~/Dropbox/QuasiAsterisk/QuasiViz/")

require("R.utils")
sourceDirectory("src/")

setwd("~/Dropbox/Projects/Martin/int/")
dir <- readDirectory("ALL_5e-4/")

pdf("entropy.pdf",15,7)
plotEntropy(dir,legendPosition="bottom")
dev.off()

pdf("hapdist.pdf",15,7)
plotHapDist(dir,ymin=5e-4)
dev.off()

mds <- computeMDS(dir)

pdf("mds.pdf",20,20)
par(mfrow=c(1,1))
plotMDSsingle(mds,7)
dev.off()

pdf("mds.pdf",30,40)
par(mfrow=c(4,3))
plotMDS(mds,offset=1,circleBorder=rgb(0,0,0,alpha=.2),cexmain=3,mar=c(0,0,3,0))
dev.off()

pdf("mds.pdf",30,40)
par(mfrow=c(4,3))
plotMDS(mds,offset=1,circleBorder=NULL,cexmain=3,mar=c(0,0,3,0),cex=1)
dev.off()

min(sapply(dir, function(x) min(x$frequencyDistribution)))

dir[[1]]$sequences