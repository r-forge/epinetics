sim.prob.epi.by.hist <-
  function(net.list, epi.cutoff=10){
  num.nets <- length(net.list)
  D <- data.frame(rep.num=1:num.nets, seed.mem=rep(99, num.nets), epi=rep(0,num.nets))
  for(i in 1:num.nets){
    net <- o2[[i]]
    patient.zero.index <- which(net %v% "time_infected" == 0)
    D$seed.mem[i] <- (net %v% "immune_memory")[patient.zero.index]
    epi.size <- sum(net %v% "status" %in% "recovered")
    if (epi.size > epi.cutoff){
      D$epi[i] <- 1
    }
  }
  D
}
