
havel.hakimi <-
  function (degree.sequence){
    num.edges <- sum(degree.sequence)/2
    heads <- vector(mode="integer", length=num.edges)
    tails <- vector(mode="integer", length=num.edges)
    x <- .C("havel_hakimi_R", as.integer(degree.sequence), 
       as.integer(length(degree.sequence)), 
       as.integer(heads), as.integer(tails), PACKAGE="epinetics")
  network(cbind(x[[3]],x[[4]]), directed=FALSE)
}

