sample.host.id.with.hist <-
  function(net, size, hist){
    hist.set <- unique(net %v% "infection_history")
    if (!is.element(hist, hist.set)){
      stop("There are no hosts with history HIST")
    }
    pool <- which(net %v% "infection_history" %in% hist)
    if(length(pool)==1)
      return(rep(pool, size))
    else{
      return(sample(pool, size=size, replace=TRUE))
    }
}
