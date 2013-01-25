readFasta <- function(input) {
  require(entropy)
  q <- scan(input, character(0))
  L <- 0
  for (i in 1:length(q)) {
    if (i %% 2 == 0) {
      L <- max(L,nchar(q[i]))
    }
  }
  m <- matrix(0,nrow=5,ncol=L)
  d <- matrix(0,nrow=length(q)/2,ncol=L)
  p <- vector(length=ceiling(length(q)/2))
  f <- vector(length=ceiling(length(q)/2))
  s <- vector(length=ceiling(length(q)/2))
  y <- vector(length=ceiling(length(q)/2))
  
  for (i in 1:length(q)) {
    if (i %% 2 == 0) {
      s <- unlist(strsplit(toupper(q[i]),""))
      for (j in 1:nchar(q[i])) {
        if (s[j] == "A") {
          d[i/2,j] <- 0
          m[1,j] <- m[1,j]+p[i/2]
        } else if (s[j] == "C") {
          d[i/2,j] <- 1
          m[2,j] <- m[2,j]+p[i/2]
        } else if (s[j] == "G") {
          d[i/2,j] <- 2
          m[3,j] <- m[3,j]+p[i/2]
        } else if (s[j] == "T") {
          d[i/2,j] <- 3
          m[4,j] <- m[4,j]+p[i/2]
        } else if (s[j] == "-") {
          d[i/2,j] <- 4
          m[5,j] <- m[5,j]+p[i/2]
        }
      }
    } else {
      header <- unlist(strsplit(q[i],"_"))
      splitLength <- length(header)
      if (splitLength >= 2) {
        p[1+(i-1)/2] <- as.double(header[2])
      }
      if (splitLength >= 3) {
        f[1+(i-1)/2] <- as.double(header[3])
      }
      if (splitLength >= 4) {
        s[1+(i-1)/2] <- as.integer(header[4])
      }
    }
  }
  
  entropyP <- apply(m,2,entropy)
  entropyP.mean <- mean(entropyP)
  entropyP.sd <- sd(entropyP)
  name <- unlist(strsplit(unlist(strsplit(input,"/"))[length(unlist(strsplit(input,"/")))],'\\.'))[1]
  list(sequences=d,
       frequencyDistribution=p,
       quasispeciesEntropy=list(data=entropyP,mean=entropyP.mean,sd=entropyP.sd),
       fitnessDistribution=f,
       significance=s,
       name=name)
}

readDirectory <- function(path) {
  require(multicore)
  files <- list.files(path)
  mclapply(files,function(x) readFasta(paste(path,x,sep="/")))
}