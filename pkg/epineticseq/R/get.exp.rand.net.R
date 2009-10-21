get.exp.rand.net <-
  function(netsize=1000, mean.degree=15){
  degs <- ceiling(rexp(netsize, rate=1/mean.degree))
  if (sum(degs) %% 2 != 0) { degs[1] <- degs[1] + 1 }
  if (max(degs) > netsize)
    stop("Degree sequence contains degree that is higher than highest possible number of neighbors")
  g <- degree.sequence.game(degs, method="vl")
  el <- get.edgelist(g) + 1
  g <- network(el, directed=FALSE)
  g
}
