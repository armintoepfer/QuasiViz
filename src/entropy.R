plotEntropy <- function(fastaList,cex = 1, main="Entropy", legendPosition = "topleft") {
  require(entropy)
  entropyList <- sapply(dir,function(x) entropy(x$frequencyDistribution))
  entropyListMean <- sapply(dir,function(x)x$quasispeciesEntropy$mean)

  cc <- cex
  l <- length(fastaList)
  par(mar=c(5.1,4.1,4.1,5.1))
  plot(0:(l-1),entropyList,pch=16,type="b",xlab="Time points",xaxt="n",ylab="Frequency distribution entropy",cex=cc+.5,cex.axis=cc,cex.lab=cc+.5,cex.main=cc+1,main=main)
  axis(1,at=0:(l-1),1:l,cex.axis=cc)
  par(new=T)
  plot(0:(l-1),entropyListMean,axes=F,pch=4,type="b",lty=2,ylab="",xlab="",cex=2,cex.axis=cc,cex.lab=cc+.5,cex.main=cc)
  axis(4,cex.axis=cc)
  mtext(side=4,text="Mean quasispecies entropy",padj=3,cex=cc+.5)
  legend(legendPosition,col=rep(1,2),lty=1:2,pch=c(16,4),c("Distribution entropy","Mean quasispecies entropy"),bty="n",cex=cc+.5)
}