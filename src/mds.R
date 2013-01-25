fastaList <- dir

computeMDS <- function(fastaList) {
  require(Rlof)
  require(mgcv)
  
  A <- matrix(0,0,1,ncol=ncol(fastaList[[1]]$sequences))
  l <- length(fastaList)
  hapSizes <- c()
  for (i in 1:l) {
    A <- rbind(A,fastaList[[i]]$sequences)
    hapSizes <- c(hapSizes,nrow(fastaList[[i]]$sequences))
  }
  Au <- uniquecombs(A)
  ind <- attr(Au,"index")
  
  d <- distmc(Au)
  fit <- cmdscale(d,2,eig=TRUE)
  
  order <- lapply(1:l,function(x) if(x==1) {1:cumsum(hapSizes)[1]} else (cumsum(hapSizes)[x-1]+1):cumsum(hapSizes)[x])
  
  xmax <- max(fit$points[,1])
  xmin <- min(fit$points[,1])
  ymax <- max(fit$points[,2])
  ymin <- min(fit$points[,2])
  
  freqMin <- min(sapply(fastaList,function(x) min(log10(x$frequencyDistribution))))
  
  points <- mclapply(1:l,function(i) list(x=fit$points[ind[unlist(order[i])],1],
                                          y=fit$points[ind[unlist(order[i])],2],
                                          freq=sapply(fastaList[[i]]$frequencyDistribution,function(x) log10(x)-freqMin+0.1),
                                          fitness=fastaList[[i]]$fitnessDistribution,
                                          fitnessSignifance=fastaList[[i]]$significance,
                                          name=fastaList[[i]]$name))
  
  list(points=points,xlim=c(xmin:xmax),ylim=c(ymin:ymax))
}

plotMDSsingle <- function(mds,index=1,cexmain=2,offset=0.5,circleBorder="gray",mar=c(0,0,4,10)) {
  require(RColorBrewer)
  require(plotrix)
  
  par(mar=mar)
  par(mar=c(0,0,4,0))
  layout(matrix(c(1,2),1,2,byrow=TRUE),width=c(1.2,.1))
  i <- index
  sizes <- mds$points[[i]]$freq
  pg <- ceiling(sizes*100)
  ccc <- rainbow(pg,alpha=.4,start = 0, end = 1)
  colors <- ccc[pg]
  symbols(mds$points[[i]]$x,mds$points[[i]]$y,
          xlim=c(min(mds$xlim)-offset,max(mds$xlim)+offset),ylim=c(min(mds$ylim)-offset,max(mds$ylim)+offset),
          circles=sizes,bg=colors,fg=circleBorder,
          main=paste("MDS",mds$points[[i]]$name),
          xlab="",ylab="",xaxt="n",yaxt="n",
          cex.main=cexmain)
  plot(0,0,type="n",axes=F,ylab="",xlab="",ylim=c(0,1),xlim=c(0,1))
  colorlegend(rainbow(pg,start = 0, end = 1),ylim=c(0,1),xlim=c(0,1),label=c(1e-4,1e-3,1e-2,1))
}

plotMDS <- function(mds,cexmain=2,offset=0.5,circleBorder="gray",mar=c(0,0,4,0)) {
  require(RColorBrewer)
  require(plotrix)
  
  par(mar=mar)
  for(i in 1:length(mds$points)) {
    sizes <- mds$points[[i]]$freq
    pg <- ceiling(sizes*100)
    ccc <- rainbow(pg,alpha=.4,start = 0, end = 1)
    colors <- ccc[pg]
    symbols(mds$points[[i]]$x,mds$points[[i]]$y,
            xlim=c(min(mds$xlim)-offset,max(mds$xlim)+offset),ylim=c(min(mds$ylim)-offset,max(mds$ylim)+offset),
            circles=sizes,bg=colors,fg=circleBorder,
            main=paste("MDS",mds$points[[i]]$name),
            xlab="",ylab="",xaxt="n",yaxt="n",
            cex.main=cexmain)
  }
}