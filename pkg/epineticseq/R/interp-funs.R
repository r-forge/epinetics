build.circulant <-
  function (N=10, mean.k=2){
    m <- matrix(0, N, N)
    K <- mean.k / 2
    if (mean.k %% 2 )    #if odd
      print(paste("warning: only even mean degrees supported"))
    for (j in 1:N){
      nb = K
      non.nb = N - 1 - 2*K
      if (j + 1 <= N){
        for ( i in (j+1):N){
          if (nb){
            m[i,j]=1
            nb = nb -1
          }
          else if (non.nb){
            m[i,j]=0
            non.nb = non.nb - 1
          }
          else{
            m[i,j]=1
          }
        }
      }
      if (j-1 >= 1)
        for ( i in 1:(j-1)){
          if (nb){
            m[i,j]=1
            nb = nb -1
          }
          else if (non.nb){
            m[i,j]=0
            non.nb = non.nb - 1
          }
          else{
            m[i,j]=1
          }
      }
   }
   if (nb || non.nb) print("error in indexing")
    else network(m, directed=FALSE)
  }
                 

`degree.pgf` <-
function (network, x, T){
    deg.set <- degree(network, cmode="indegree")
    deg.dist <- table(deg.set)
    deg.dist.norm <- deg.dist/sum(deg.dist)
    sum <- 0
    for(i in 1:length(deg.dist.norm)){
      k <- as.numeric(names(deg.dist.norm)[i])
      pk <- as.numeric(deg.dist.norm[i])
      sum <- sum + pk*(1+(x-1)*T)^(k)
    }
    sum
    
  }

`excess.degree.pgf` <-
function (network, x, T){
    deg.set <- degree(network, cmode="indegree")
    mean.k <- mean(deg.set)
    deg.dist <- table(deg.set)
    deg.dist.norm <- deg.dist/sum(deg.dist)
    sum <- 0
    for(i in 1:length(deg.dist.norm)){
      k <- as.numeric(names(deg.dist.norm)[i])
      pk <- as.numeric(deg.dist.norm[i])
      sum <- sum + k*pk*(1+(x-1)*T)^(k-1)
    }
    sum/mean.k
    
  }

`final.epi.size` <-
function(network, T, tol=0.001, max.iterations=10){
    u <- 0
    diff <- 1
    counter <- 0
    while(diff>tol && counter < max.iterations){
      u.new <- excess.degree.pgf(network=network, T=T, x=u)
      diff <- abs(u.new-u)
      u <- u.new
      counter <- counter + 1
    }
    
    network.size(network)*(1-degree.pgf(network=network, x=u, T=T))
  }

`get.mean.outbreak` <-
function (network, T){
    if (is.directed(network))
      print("warning: this formula is only for undirected networks")
    k <- mean(degree(network, cmode="indegree"))
    k.sqd <- mean(degree(network, cmode="indegree")^2)
    G.1.prime <- (k.sqd - k)/k
    G.0.prime <- k
    mos <- 1 + T*G.0.prime/(1-T*G.1.prime)
    mos
  }

`get.t.crit` <-
function (network){
    if (is.directed(network))
      print("warning: this formula is only for undirected networks")
    k <- mean(degree(network, cmode="indegree"))
    k.sqd <- mean(degree(network, cmode="indegree")^2)
    G.1.prime <- (k.sqd - k)/k
    t.crit <- 1/G.1.prime
    t.crit
  }

`immunize.deg.dist.tail` <-
function(numhosts, strain.id, upper=TRUE){
    degdist <- degree(epi.network)
    if (upper) lapply(c(as.numeric(names(degdist[order(degdist, decreasing=TRUE)][1:numhosts]))), immunize, strain.id)
    else lapply(c(as.numeric(names(degdist[order(degdist, decreasing=FALSE)][1:5]))), immunize, strain.id)
  }

`immunize` <-
function (ego.id, strain.id){
    set.vertex.attribute (epi.network, "immune_memory", strain.id, ego.id)
  }

`infect` <-
function (nb.id, strain.id){
    set.vertex.attribute (epi.network, "status", "infectious", nb.id)
    set.vertex.attribute (epi.network, "time_infected", sim.time, nb.id)
    set.vertex.attribute (epi.network, "infection_history", strain.id, nb.id)
    infectious.ids <<- append (infectious.ids, nb.id)
#      print (paste ("total.rates before infection of", nb.id, "was", round(total.rates,3)))    
    hood <- get.neighborhood (epi.network, nb.id)
    nb2.states <- get.vertex.attribute (epi.network, "status")[hood]
    sus.nb2.ids <- hood[ which ( nb2.states %in% "susceptible" )]
    #add event of strain mutating
    cl <- call ("strain.mutate", nb.id)
    key <- paste (nb.id, "strain.mutates", sep="")
    total.rates <<- total.rates + mutation.rate
    event.list[[key]] <<- c(cl, mutation.rate)
    #add event of recovery
    cl <- call ("recovery", nb.id)
    key <- paste (nb.id, "recovers", sep="")
    total.rates <<- total.rates + recovery.rate
    event.list[[key]] <<- c(cl, recovery.rate)
    for (y in seq (along=sus.nb2.ids)){  #add events of the newly infectious host infecting his S nbs
      nb2.id <- sus.nb2.ids[y]
      cl <- call ("infect", nb2.id, strain.id)
      key <- paste (nb.id, "infects", nb2.id, "with", strain.id, sep="")

      #infection rate multiplies by (1- \tau_distance), where \tau_distance = exp(-a*distance)
      nb2.immune.memory <- get.vertex.attribute (epi.network, "immune_memory")[nb2.id]
      distance <- abs (strain.id - nb2.immune.memory)
      tau <- exp(-a*distance)
      rate <- infection.rate * (1 - tau)
      total.rates <<- total.rates + rate 
      event.list[[key]] <<- c(cl, rate)
#      print (paste (key, "distance:", distance,  " rate:", rate))
    }

    #remove events of the newly infected host getting infected by any other infectious nb
    inf.nb2.ids <- hood[ which ( nb2.states %in% "infectious" )]
    for (y in seq (along=inf.nb2.ids)){
      nb2.id <- inf.nb2.ids[y]
      nb2.strain.id <- get.vertex.attribute (epi.network, "infection_history")[nb2.id]
      key <- paste (nb2.id, "infects", nb.id, "with", nb2.strain.id, sep="")
      total.rates <<- total.rates - event.list[[key]][[2]]
      event.list[[key]] <<- NULL
    }
#    print (paste ("total.rates after infection of", nb.id, "was", round(total.rates,3)))

  }

`make.multi.time.series` <-
function (sim.result, time.steps=100){
    max.strain.id <- max(na.omit(sim.result %v% "infection_history"))
    min.strain.id <- min(na.omit(sim.result %v% "infection_history"))
    mts <- make.time.series(sim.result, all.strains=FALSE, strain.vector=min.strain.id, time.steps=time.steps)
    vec <- (min.strain.id+1):max.strain.id
    for (id in vec){
      next.ts <- make.time.series(sim.result, all.strains=FALSE, strain.vector=id, time.steps=time.steps) 
      mts <- cbind(mts, next.ts[,2])
    }
    colnames(mts) <- c("time", paste("I", min.strain.id:max.strain.id, sep="."))
    mts
  }

`make.time.series` <-
function (time.steps=100, sim.result, all.strains=TRUE, strain.vector=NA){
    end.time <- max (na.omit (sim.result %v% "time_recovered"))
    start.time <- min (na.omit (sim.result %v% "time_infected"))
    step.size <- (end.time - start.time)/time.steps
    time.sequence <- seq(start.time, end.time, by=step.size)
    df <- data.frame(time=time.sequence)
    if (!all.strains){
      subset.ids <- which(sim.result %v% "infection_history" %in% strain.vector==TRUE)
      for (x in seq(along=time.sequence)){
        t <- time.sequence[x]
#       print(paste("time:",t))
        nisf <- sum((sim.result %v% "time_infected" < t)[subset.ids])
#       print(paste("number infected so far:", nisf))
        nrsf <- sum((sim.result %v% "time_recovered" < t)[subset.ids])
#       print(paste("number recovered so far:", nrsf))
        nsi <- nisf-nrsf
        df[x,2] <- nsi
      }
    }
    else{
      for (x in seq(along=time.sequence)){
        t <- time.sequence[x]
#       print(paste("time:",t))
        nisf <- sum(na.omit(sim.result %v% "time_infected" < t))
#       print(paste("number infected so far:", nisf))
        nrsf <- sum(na.omit(sim.result %v% "time_recovered" < t))
#       print(paste("number recovered so far:", nrsf))
        nsi <- nisf-nrsf
        df[x,2] <- nsi
      }
    }
    if (is.na(strain.vector) || all.strains!=TRUE){
      colnames(df)[2] <- paste("I",strain.vector, sep=".")
    }
    else{
      colnames(df)[2] <- "infectious"
    }
    df
}

`mean.epidemics.series` <-
function (network, num.reps=10, num.points=10, make.plot=0, max.trans=0.6){
    t.crit <- get.t.crit(network=network)
    if(make.plot!=0)
      curve(final.epi.size(network=network, T=x), from=t.crit, to=1)
    df <- data.frame()
    for(point in 1:num.points){
      print(paste("point",point,"of",num.points))
      trans <- ((point-1)/num.points)*(max.trans-t.crit)+t.crit - 0.01
      sizes <- list()
      for(rep in 1:num.reps){
        sim.res <- simulate.epidemic(network=network, initial.infectious=c(1:5), T=trans)
        sizes[[rep]] <- sum(sim.res %v% "status" %in% "recovered")
      }
     df <- rbind(df, c(trans, mean(as.numeric(sizes))))
    }
    if(make.plot!=0)
      points(df)
    df
  }

`mean.outbreaks.series` <-
function (network, num.reps=10, num.points=10, make.plot=0){
    t.crit <- get.t.crit(network=network)
    if(make.plot!=0)
      curve(get.mean.outbreak(network=network, T=x), from=0, to=c(t.crit-t.crit/10))
    df <- data.frame()
    for(point in 1:num.points){
      print(paste("point",point,"of",num.points))
      trans <- (point/num.points)*(9*t.crit/10)
      sizes <- list()
      for(rep in 1:num.reps){
        sim.res <- simulate.epidemic(network=network, initial.infectious=c(1), T=trans)
        sizes[[rep]] <- sum(sim.res %v% "status" %in% "recovered")
      }
     df <- rbind(df, c(trans, mean(as.numeric(sizes))))
     if(make.plot!=0)
       points(trans, mean(as.numeric(sizes)))
    }
    df
  }

`mts.stats` <-
function(mts){
    stats <- data.frame(time=mts[,1])
    counts <- list()
    for(i in 1:dim(mts)[1]){
      counts[[i]] <- vector(mode="numeric",length=0)
      for(j in 2:dim(mts)[2]){
        counts[[i]] <- append(counts[[i]], rep(j-2,mts[i,j]))
      }
    }
     
    for (i in 1:dim(mts)[1]){
      stats[i,2] <- var(counts[[i]])
      stats[i,3] <- mean(counts[[i]])
    } 
    #consitency check
    ms <- numeric(0)
    for (i in 1:dim(mts)[1]) stats[i,4] <- weighted.mean(w=mts[i,-1],x=0:(dim(mts)[2]-2))
#    if (sum( abs(stats[,3]-stats[,4]))!=0) print("means are inconsistent")
    colnames(stats)[-1] <- c("variance", "mean", "mean2")
    stats
  }

`new.epi.env` <-
function (network, T=0.3, infectious.period=5, event.list=list(), sim.time=0, infectious.ids=integer(0), status="susceptible", time.infected=as.numeric(NA), total.rates=0, mu=0.05){
    network %n% "transmissability" <- T
    infection.rate <<- -T / (infectious.period * (T -1 ))
#    infection.rate <<- -log(1-T)*(1/infectious.period)
    recovery.rate <<- 1/infectious.period
    mutation.rate <<- mu
    network %n% "infectious_period" <- infectious.period
    event.list <<- event.list
                                        #    network %n% "simulation_time" <- sim.time
    sim.time <<- sim.time
    a <<- 0.1
    infectious.ids <<- infectious.ids
    network %v% "status" <- status
    network %v% "time_infected" <- time.infected
    network %v% "immune_memory" <- -100
    network %v% "time_recovered"<- as.numeric(NA)
    network %v% "infection_history" <- integer(0)
    total.rates <<- total.rates    
    epi.network <<- network
  }

`recovery` <-
function (ego.id){
    inf.time <- get.vertex.attribute (epi.network, "time_infected")[ego.id]
    strain.id <- get.vertex.attribute (epi.network, "infection_history")[ego.id]
#    print (paste ("host", ego.id, "infected at:", round(inf.time,2), "is recovering at", round(sim.time, 2), "the total.rates is:", round(total.rates,2)))
    set.vertex.attribute (epi.network, "status", "recovered", ego.id)
    set.vertex.attribute (epi.network, "time_recovered", sim.time, ego.id)
    infectious.ids <<- infectious.ids[-which(infectious.ids==ego.id)]
    #remove event of this host infection mutation to a new variant
    key <- paste (ego.id, "strain.mutates", sep="")
    total.rates <<- total.rates - event.list[[key]][[2]]
    event.list[[key]] <<- NULL
    #remove event of this host recovering from event.list
    key <- paste (ego.id, "recovers", sep="")
    total.rates <<- total.rates - event.list[[key]][[2]]
    event.list[[key]]<<- NULL

    hood <- get.neighborhood (epi.network, ego.id)
    nb.states <- get.vertex.attribute (epi.network, "status")[hood]
    sus.nb.ids <- hood[ which ( nb.states %in% "susceptible" )]
#    print ( paste("susceptible nbs are", sus.nb.ids))
    for (y in seq (along=sus.nb.ids)){  #remove events of the newly recovered host infecting his S nbs
      
      nb.id <- sus.nb.ids[y]
      key <- paste (ego.id, "infects", nb.id, "with", strain.id, sep="")
#      print (paste("key being removed is", key))
      total.rates <<- total.rates - event.list[[key]][[2]]
      event.list[[key]] <<- NULL
    }
#    print (paste ("total.rates now is:", round(total.rates,2)))
  }

`simulate.epidemic` <-
function (network, max.steps=100000, verbose=FALSE, tolerance=0.0001, initial.infectious=c(6,2,3,4,5), pre.immunize=FALSE, initial.immune=c(6,7), tail=FALSE, upper.tail=TRUE, T=0.3){
    new.epi.env (network, infectious.period=10, T=T)
    if (tail) immunize.deg.dist.tail (3, -1, upper=upper.tail)
    else if (pre.immunize) lapply (initial.immune, immunize, -1)   
    lapply (initial.infectious, infect, 0)
    steps <- 0
#    initializations ()
    if (verbose){
      print ("step | time | number infectious")
      while (steps <= max.steps && total.rates > tolerance){
        steps <- steps + 1
        simulate.next.event()
        print (paste (steps, round(sim.time,3), length(infectious.ids), round(total.rates,3), sep=" | "))
      }
    }
    else{
      while (steps <= max.steps && total.rates > tolerance){
      steps <- steps + 1
      simulate.next.event()
      }
    }
    epi.network
  }

`simulate.next.event` <-
function (){
    rand.one <- runif(1)
    time.interval <- -log(rand.one)/total.rates
    sim.time <<- sim.time + time.interval
    rand.two <- runif(1, 0, total.rates)
    rate.sum <- 0
    for (x in seq (along=event.list)){
      rate.sum <- rate.sum + event.list[[x]][[2]]
      if (rate.sum >= rand.two){
        break
      }
    }
    eval(event.list[[x]][[1]])

  }

ssa.epidemics.series <-
  function (network, num.reps = 10, num.points = 10, make.plot = 0, 
            max.trans = 0.6) 
{
  t.crit <- get.t.crit(network = network)
  if (make.plot != 0) 
    curve(final.epi.size(network = network, T = x), from = t.crit, 
          to = 1)
  df <- data.frame()
  for (point in 1:num.points) {
    print(paste("point", point, "of", num.points))
    trans <- ((point - 1)/num.points) * (max.trans - t.crit) + 
      t.crit - 0.01
    sizes <- list()
    for (rep in 1:num.reps) {
      sim.res <- epi.sim.ssa(net = network, trans = trans, max.steps=100000)
      sizes[[rep]] <- sum(sim.res %v% "status" %in% "recovered")
    }
    df <- rbind(df, c(trans, mean(as.numeric(sizes[sizes>100]))))
  }
  if (make.plot != 0) 
    points(df)
  df
}
ssa.outbreaks.series <-
  function (network, num.reps = 10, num.points = 10, make.plot = 0) 
{
  t.crit <- get.t.crit(network = network)
  if (make.plot != 0) 
    curve(get.mean.outbreak(network = network, T = x), from = 0, 
          to = c(t.crit - t.crit/10))
  df <- data.frame()
  for (point in 1:num.points) {
    print(paste("point", point, "of", num.points))
    trans <- (point/num.points) * (9 * t.crit/10)
    sizes <- list()
    for (rep in 1:num.reps) {
      sim.res <- epi.sim.ssa(net = network, trans = trans, max.steps=100000)
      sizes[[rep]] <- sum(sim.res %v% "status" %in% "recovered")
    }
    df <- rbind(df, c(trans, mean(as.numeric(sizes))))
    if (make.plot != 0) 
      points(trans, mean(as.numeric(sizes)))
  }
  df
}
`strain.mutate` <-
function (ego.id){
#    print ("mutations")
    old.strain <- get.vertex.attribute (epi.network, "infection_history")[ego.id]
    new.strain <- old.strain + 1
    set.vertex.attribute (epi.network, "infection_history", new.strain, ego.id)

    hood <- get.neighborhood (epi.network, ego.id)
    nb.states <- get.vertex.attribute (epi.network, "status")[hood]
    sus.nb.ids <- hood[ which ( nb.states %in% "susceptible" )]
    for (y in seq (along=sus.nb.ids)){  #update the infection calls 
      nb.id <- sus.nb.ids[y]
      keyold <- paste (ego.id, "infects", nb.id, "with", old.strain, sep="")
#      print (paste("key being removed is", keyold, "total.rates:", total.rates))
      total.rates <<- total.rates - event.list[[keyold]][[2]]
      event.list[[keyold]] <<- NULL
      
      keynew <- paste (ego.id, "infects", nb.id, "with", new.strain, sep="")
#      print (paste("key being added is", keynew, "total.rates", total.rates))
      cl <- call ("infect", nb.id, new.strain)

      #infection rate multiplies by (1- \tau_distance), where \tau_distance = exp(-a*distance)
      nb.immune.memory <- get.vertex.attribute (epi.network, "immune_memory")[nb.id]
      distance <- abs (new.strain - nb.immune.memory)
      tau <- exp(-a*distance)
      rate <- infection.rate * (1 - tau)
      total.rates <<- total.rates + rate 
      event.list[[keynew]] <<- c(cl, rate)
#     print (paste (keynew, "distance:", distance,  " rate:", rate))
    }  
  }

`ts.stats.plotter` <-
function (mts, mu=0.05){
    stats2 <- mts.stats(mts)
    plot(stats2[,1],stats2[,1]*mu, xlab="days", ylab="")
    points(stats2[-c(1,101),1],stats2[-c(1,101),2], pch="v")
    points(stats2[-c(1,101),1],stats2[-c(1,101),3], pch="m")
  }

vert.status.plot <-
   function (net){
   plot(net, displaylabels=TRUE, vertex.col=match(a%v% "status", levels(as.factor((a %v% "status"))) ))
}
