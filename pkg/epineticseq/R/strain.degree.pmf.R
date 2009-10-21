strain.degree.pmf <-
  function (net) 
{
  types <- unique(net %v% "infection_history")
  coded.history <- (net %v% "infection_history") + 1
  max.types <- max(na.omit(coded.history)) + 1
  coded.history[is.na(coded.history)] <- max.types
  type.codes <- unique(coded.history)
#  print(paste("types:", types))
#  print(paste("type.codes:", type.codes))
  pmf.list <- list()
#  for (i in 1:1){
  for (i in 1:length(types)) {
    cur.code <- type.codes[i]
    set <- which(coded.history == cur.code)
#    print(paste("set:", set))
    l <- lapply(set, function(x) net$val[[x]]$strain_degrees)
#    print(l)
    l.t <- list()
    for (j in 1:max.types) {
#      print(j)
      l.t[[paste("sc", j, sep = "")]] <- rapply(lapply(l, 
                                                       function(x) x[[j]]), function(x) x)
    }
#    print(l.t)
    tab <- table(l.t)
    tab <- tab/sum(tab)
    pmf.list[[cur.code]] <- tab
  }
  pmf.list
}
