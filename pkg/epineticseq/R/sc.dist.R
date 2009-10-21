sc.dist <-
  function(net){
    codes <- net %v% "infection_history"
    codes <- codes + 1
    codes[is.na(codes)] <- max(codes[!is.na(codes)])+1
    x <- tabulate(codes )
    x <- x/sum(x)
    x
  }
