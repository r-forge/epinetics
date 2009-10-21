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

