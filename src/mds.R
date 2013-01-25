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

plotMDS <- function(mds,cexmain=2,offset=0.5,circleBorder="gray") {
  require(RColorBrewer)
  require(plotrix)
  pm_old <- par()$mar
  
  par(mar=c(0,0,4,0))
  for(i in 1:length(mds$points)) {
    sizes <- mds$points[[i]]$freq
    pg_max <- ceiling(max(sizes)*100)
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