make.time.series <-
  function (sim.result, time.steps = 100, all.strains = TRUE, strain.vector = NA) 
{
  end.time <- max(na.omit(sim.result %v% "time_recovered"))
  start.time <- min(na.omit(sim.result %v% "time_infected"))
  step.size <- (end.time - start.time)/time.steps
  time.sequence <- seq(start.time, end.time, by = step.size)
  df <- data.frame(time = time.sequence)
  if (!all.strains) {
    subset.ids <- which(sim.result %v% "infection_history" %in% 
                        strain.vector == TRUE)
    for (x in seq(along = time.sequence)) {
      t <- time.sequence[x]
      nisf <- sum((sim.result %v% "time_infected" < t)[subset.ids])
      nrsf <- sum((sim.result %v% "time_recovered" < t)[subset.ids])
      nsi <- nisf - nrsf
      df[x, 2] <- nsi
    }
  }
  else {
    for (x in seq(along = time.sequence)) {
      t <- time.sequence[x]
      nisf <- sum(na.omit(sim.result %v% "time_infected" < 
                          t))
      nrsf <- sum(na.omit(sim.result %v% "time_recovered" < 
                          t))
      nsi <- nisf - nrsf
      df[x, 2] <- nsi
    }
  }
  if (is.na(strain.vector) || all.strains != TRUE) {
    colnames(df)[2] <- paste("I", strain.vector, sep = ".")
  }
  else {
    colnames(df)[2] <- "infectious"
  }
  df
}
