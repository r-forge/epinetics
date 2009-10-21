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

