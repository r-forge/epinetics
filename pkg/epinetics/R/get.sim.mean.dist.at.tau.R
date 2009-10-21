`get.sim.mean.dist.at.tau` <-
function(net, tau){
    ts <- net %n% "time_series"
    if(length(ts)==1)
      stop("Network does not have full time series")
    deltas <- abs(ts - tau)
    index <- which(deltas==min(deltas))
    closest.time <- ts[index]
    mean.dist <- (net %n% "mean_distance_series")[index]
    cbind(closest.time, mean.dist)
  }

