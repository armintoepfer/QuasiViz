setwd("~/Dropbox/QuasiAsterisk/QuasiViz/")

require("R.utils")
sourceDirectory("src/")

#readFasta("~/Dropbox/Projects/Martin/int/20130117/ALL/08-16621.fasta")
dir <- readDirectory("data/")

plotEntropy(dir,legendPosition="bottom")
plotHapDist(dir)

mds <- computeMDS(dir)

pdf("mds.pdf",20,20)
par(mfrow=c(1,1))
plotMDS(mds)
dev.off()

pdf("mds.pdf",40,30)
par(mfrow=c(3,4))
plotMDS(mds,offset=1,circleBorder="black")
dev.off()
