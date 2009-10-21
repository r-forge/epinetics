make.var.in.mean.vs.t.df <-
  function(network, reps=20, trannies=c(0.09, 0.17, 0.26, 0.34), nu=0.2, kappa=14*.96,...){
    df <- data.frame(T=trannies,
                     sigma.tau=array(-99, dim=length(trannies)),
                     sigma.inf=array(-99, dim=length(trannies)))
    for( i in 1:length(trannies)){
      T <- trannies[i]
      o <- epi.sim.ssa(network, sample.initial.infectious=TRUE,
                       trans=T, reps=reps, infectious.period=1/nu,
                       make.time.series=1,...)
      df$sigma.inf[i] <- var(as.numeric(lapply(o, function(x) mean(na.omit(x %v% "infection_history")))))
      N <- network.size(network)
      tau <- get.tau(trans=T, nu=nu, kappa=kappa, N=N)
      temp <- mean.dist.tau.df(o, tau)
      df$sigma.tau[i] <- var(temp$mdists)
    }
    df
}
