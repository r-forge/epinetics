`mts.stats` <-
function(mts){
    stats <- data.frame(time=mts[,1])
    counts <- list()
    for(i in 1:dim(mts)[1]){
      counts[[i]] <- vector(mode="numeric",length=0)
      for(j in 2:dim(mts)[2]){
        counts[[i]] <- append(counts[[i]], rep(j-2,mts[i,j]))
      }
    }
     
    for (i in 1:dim(mts)[1]){
      stats[i,2] <- var(counts[[i]])
      stats[i,3] <- mean(counts[[i]])
    } 
    #consitency check
    ms <- numeric(0)
    for (i in 1:dim(mts)[1]) stats[i,4] <- weighted.mean(w=mts[i,-1],x=0:(dim(mts)[2]-2))
#    if (sum( abs(stats[,3]-stats[,4]))!=0) print("means are inconsistent")
    colnames(stats)[-1] <- c("variance", "mean", "mean2")
    stats
  }

