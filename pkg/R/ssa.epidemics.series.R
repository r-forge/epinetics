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
