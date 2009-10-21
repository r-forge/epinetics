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

