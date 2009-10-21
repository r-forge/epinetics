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

