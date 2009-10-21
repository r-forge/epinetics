`immunize.deg.dist.tail` <-
function(numhosts, strain.id, upper=TRUE){
    degdist <- degree(epi.network)
    if (upper) lapply(c(as.numeric(names(degdist[order(degdist, decreasing=TRUE)][1:numhosts]))), immunize, strain.id)
    else lapply(c(as.numeric(names(degdist[order(degdist, decreasing=FALSE)][1:5]))), immunize, strain.id)
  }

