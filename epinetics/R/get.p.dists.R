get.dist.vs.time.p <-
  function(net, sample.frac=0.05){
    infected.indices <- which(net %v% "status" %in% "recovered")
    total.infected <- length(infected.indices)
    if(total.infected <= 0){
      stop("No recovered individuals in this network\n")
    }
    sample.size <- round(total.infected*sample.frac, digits=0)
    sample.indices <- infected.indices[sample(x=total.infected, size=sample.size, replace=FALSE)]
    distances <- (net %v% "infection_history")[sample.indices]
    recovery.times <- (net %v% "time_recovered")[sample.indices]
    mod <- lm(distances ~ 0+recovery.times)
    cf <- coef(summary(mod))
    p <- cf["recovery.times","Pr(>|t|)"]
    p
  }
    
    
    
