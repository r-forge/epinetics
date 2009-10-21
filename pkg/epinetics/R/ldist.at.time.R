ldist.at.time <-
  function(net.list, time=9.7){
    rapply(lapply(net.list, dist.at.time, time=time), function(x) x)
  }
