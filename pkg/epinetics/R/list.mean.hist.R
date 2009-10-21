list.mean.hist <-
  function(net.list){
    rapply(lapply(net.list, function(x) mean(na.omit(x %v% "infection_history"))), function(x) x)
  }
