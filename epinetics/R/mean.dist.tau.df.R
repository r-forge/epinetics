mean.dist.tau.df <-
  function(net.list, tau){
    l <- lapply(net.list, get.sim.mean.dist.at.tau, tau)
    data.frame(times=rapply(l, function(x) x[1]),
               mdists=rapply(l, function(x) x[2]))
  }
    
