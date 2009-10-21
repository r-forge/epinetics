get.infection.tree <- function(network){
  infection.edges <- data.frame(infector.id=sapply(network$val, function(x) x$infector_id),
                                ego.id=sapply(network$val, function(x) x$vertex.names))
  patient.zero <- which(infection.edges$infector.id==0)
  infection.edges$infector.id[patient.zero] <- infection.edges$ego.id[patient.zero]
  network(na.omit(infection.edges))
}
