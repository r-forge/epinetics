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

