ll.crunch <-
function(list.of.list.of.nets, function.name){
  ll <- list.of.list.of.nets
  lapply(ll, function(x) {
    if(is.list(x) && !is.network(x))
      lapply(x, function.name)
  })
}
