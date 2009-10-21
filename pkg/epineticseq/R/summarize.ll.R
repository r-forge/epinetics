summarize.ll <- function(list.of.list.of.nets){
  ll <- list.of.list.of.nets
  lapply(ll, function(x) {
    if(is.list(x) && !is.network(x))
      lapply(x, function(x) sum("recovered" == x %v% "status"))
  })
}
       
