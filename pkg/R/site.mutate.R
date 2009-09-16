site.mutate <-
  function(sequence.id, position, base){
    sequence.list[[sequence.id]][position] <<- base
  }
