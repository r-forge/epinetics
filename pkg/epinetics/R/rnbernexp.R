rnbernexp <- function(n, nv, p = 0.5, onset.hazard = 1, termination.hazard = 1){
    nets <- list()
    for(i in 1:n)
      nets[[i]] <- .Call("rnbernexp_R", network.initialize(nv, directed = FALSE), p, onset.hazard, termination.hazard, PACKAGE = "epinetics")
    if(i > 1)
      nets
    else
      nets[[1]]
  }

