`get.mean.outbreak` <-
function (network, T){
    if (is.directed(network))
      print("warning: this formula is only for undirected networks")
    k <- mean(degree(network, cmode="indegree"))
    k.sqd <- mean(degree(network, cmode="indegree")^2)
    G.1.prime <- (k.sqd - k)/k
    G.0.prime <- k
    mos <- 1 + T*G.0.prime/(1-T*G.1.prime)
    mos
  }

