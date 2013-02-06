plotHapDist <- function(fastaList,cex = 3, lwd=3, ymin=1e-3) {
  par(mar=c(5,6,4,1))
  l <- length(fastaList)
  xmax <- max(sapply(fastaList,function(x) length(x$frequencyDistribution)))
  
  plot(1,1,type="n",log="yx",cex.axis=cex,xaxt="n",cex.main=2,axes=F,xlim=c(1,xmax),yaxt="n",ylim=c(ymin,1),xlab="",ylab="",main="Haplotype distributions")
  mtext(side=2,text="Frequency",cex=2,padj=-3.5)
  mtext(side=1,text="Haplotypes",cex=2,padj=3)
  axis(2, at=c(1e-3,1e-2,1e-1,1), labels=c(expression(10^-3),expression(10^-2),expression(10^-1),1), lwd.ticks=1,cex.axis=1.5)
  axis(1)
  
  for (i in 1:l) {
    x <- fastaList[[i]]$frequencyDistribution
    lines(1:length(x),x,lwd=lwd,col=rainbow(l)[i])
  }
  legend("topright",bty="n",col=rainbow(l),lwd=lwd,lty=rep(1,l),legend=sapply(fastaList,function(x) x$name))
}



#axis(2, at=c(1e-3,1e-2,1e-1,1), labels=c(expression(10^-3),expression(10^-2),expression(10^-1),1), lwd.ticks=1,cex.axis=1.5)
#legend("topright",bty="n",col=rainbow(11),lwd=3,lty=rep(1,11),legend=c("Time point 1","Time point 2","Time point 3","Time point 4","Time point 5","Time point 6","Time point 7","Time point 8","Time point 9","Time point 10","Time point 11"))
