list.mean.dist <-
  function(net.list, tau){
    data.frame(time=rapply(lapply(o, get.sim.mean.dist.at.tau, tau=tau), function(x) x[1]),  mean.dist=rapply(lapply(o, get.sim.mean.dist.at.tau, tau=tau), function(x) x[2]))
  }
