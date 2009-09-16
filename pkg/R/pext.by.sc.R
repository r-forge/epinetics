pext.by.sc <-
  function(net, reps=10, trans=0.7, epi.cutoff=10,...){
    s.set <- unique(net %v% "infection_history")
    E <- data.frame(last.infection=s.set, pext=0)
    for (i in 1:length(s.set)){
      s <- s.set[i]
      ids <- sample.host.id.with.hist(net, size=reps, hist=s)
      o2 <- (lapply(ids, function(x) epi.sim.ssa(net=net, reps=1, trans=trans, initial.infectious.id = x, mutation.rate=0,...)))
      a <- lapply(o2, function(x) fs <- sum(x %v% "status" %in% "recovered"))
      a <- lapply(a, function(x) {if(x > epi.cutoff) 1 else 0})
      b <- 1 - mean(as.numeric(a))
      E[i,2] <- b
    }
    E
  }
