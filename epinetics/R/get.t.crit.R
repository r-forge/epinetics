`get.t.crit` <-
function (network){
    if (is.directed(network))
      print("warning: this formula is only for undirected networks")
    k <- mean(degree(network, cmode="indegree"))
    k.sqd <- mean(degree(network, cmode="indegree")^2)
    G.1.prime <- (k.sqd - k)/k
    t.crit <- 1/G.1.prime
    t.crit
  }

