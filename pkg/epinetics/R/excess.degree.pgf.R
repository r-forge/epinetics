`excess.degree.pgf` <-
function (network, x, T){
    deg.set <- degree(network, cmode="indegree")
    mean.k <- mean(deg.set)
    deg.dist <- table(deg.set)
    deg.dist.norm <- deg.dist/sum(deg.dist)
    sum <- 0
    for(i in 1:length(deg.dist.norm)){
      k <- as.numeric(names(deg.dist.norm)[i])
      pk <- as.numeric(deg.dist.norm[i])
      sum <- sum + k*pk*(1+(x-1)*T)^(k-1)
    }
    sum/mean.k
    
  }

