add.mutation.events <-
  function(base, position, rate.mat, sequence.id){
    row.num <- match(base, rownames(rate.mat))
    row <- rate.mat[row.num,]
    for(i in 1:4){
      rate <- row[i]
      new.base <- names(row[i])
      cl <- call("site.mutate", sequence.id, position, new.base)
      key <- paste(base, position, new.base, sep="")
      total.rates <<- total.rates + rate
      event.list[[key]] <<- c(cl, rate)
    }
  
}

