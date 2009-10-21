get.reg.rand.graph <-
  function(mean.degree=15, num.nodes=1000){
    require(igraph)
    g <- degree.sequence.game(rep(mean.degree, num.nodes), method="vl")
    el <- get.edgelist(g)+1
    g <- network(el, directed=FALSE)
    g
  }
