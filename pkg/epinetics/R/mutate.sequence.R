mutate.sequence <-
function(node.id, rate.mat){
    s <- sequence.list[[node.id]]
    position <<- 1
    event.list <<- list()
    total.rates <<- 0
    sapply(s, function(x) {
      add.mutation.events(x, position, rate.mat, node.id)
      position <<- position + 1
    })
    rand <- runif(1, 0, total.rates)
    rate.sum <- 0
    for (x in seq (along=event.list)){
      rate.sum <- rate.sum + event.list[[x]][[2]]
      if (rate.sum >= rand){
        break
      }
    }
    eval(event.list[[x]][[1]])
  }

