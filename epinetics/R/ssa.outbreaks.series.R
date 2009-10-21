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
