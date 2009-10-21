get.phylo.network <- function(epi.model){
#  t <- get.infection.tree(epi.model)
  n.less.one <- max(na.omit(o %v% "phylo_id"))-1
  phylo.list <- lapply(epi.model$val, function(x) x$phylo_id)
  phylo.edge.list <- matrix(0, nrow=n.less.one, ncol=2)
  i <- 1
  lapply(phylo.list, function(x){
    j <- 1
    while(j < length(x)){
      phylo.edge.list[i,] <<- c(x[j], x[j+1])
      i <<- i+1
      j <- j+1
    }
  })
  network(phylo.edge.list)
}
