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

