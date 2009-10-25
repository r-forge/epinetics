within.epi.sim.ssa <-
  function (net, reps = 1, trans = 0.5, max.steps = 1e+05, infectious.period = 5, 
            mutation.rate = 0.05, verbose = 0, immune.decay = 0.1, initial.susceptibility = 0.99999, 
            make.time.series = 0, sample.inoc = FALSE, inoc.strain = 0, 
            sample.initial.infectious = FALSE, initial.infectious.id = 1, inoc.size=100, eff=.1, conv=.1) 
{
  on.exit(.C("simulator_free"))
  result.nets <- list()
  if (!sample.initial.infectious) {
    num.nodes <- network.size(net)
    if (initial.infectious.id > num.nodes || initial.infectious.id < 
        1) {
      stop("INITIAL.INFECTIOUS.ID must be in [1, NETWORK_SIZE].")
    }
  }
  if (initial.susceptibility == 1) {
    initial.distance = -log(1 - 0.99999)/immune.decay
    print("Warning: Approximating initial susceptibility = 1 with 0.99999")
  }
  else if (initial.susceptibility > 1 || initial.susceptibility < 
           0) {
    print("Warning: Illegal initial susceptibility, using 0.99999 instead")
    initial.distance = -log(1 - 0.99999)/immune.decay
  }
  else initial.distance <- -log(1 - initial.susceptibility)/immune.decay
  for (i in 1:reps) {
    g <- net
    if (is.na(match("recovered", g %v% "status"))) {
      g %v% "immune_memory" <- as.integer(-initial.distance)
      g %n% "this_years_strain_0_was_last_years" <- as.integer(initial.distance)
    }
    else if (sample.inoc) {
      temp <- g %v% "infection_history"
      pool <- as.integer(na.omit(temp))
      inoc <- sample(pool, 1, replace = TRUE)
      g %v% "last_years_immune_memory" <- as.integer(g %v% 
                                                     "immune_memory")
      temp[is.na(temp)] <- -inoc + (g %v% "last_years_immune_memory")[is.na(temp)]
      g %v% "immune_memory" <- as.integer(-abs(temp - inoc))
      g %n% "this_years_strain_0_was_last_years" <- inoc
    }
    else {
      if (inoc.strain < 0) {
        stop("Cannot innoculate with a negative strain")
      }
      inoc <- inoc.strain
      temp <- g %v% "infection_history"
      g %v% "last_years_immune_memory" <- as.integer(g %v% 
                                                     "immune_memory")
      temp[is.na(temp)] <- -inoc + (g %v% "last_years_immune_memory")[is.na(temp)]
      g %v% "immune_memory" <- as.integer(-abs(temp - inoc))
      g %n% "this_years_strain_0_was_last_years" <- inoc
    }
    g %n% "replicate_number" <- i
    g %n% "time_series" <- as.double(0)
    g %n% "case_series" <- as.integer(1)
    g %n% "cov_series" <- as.double(0)
    g %n% "fitness_series" <- as.double(0)
    g %n% "mean_distance_series" <- as.double(0)
    g %v% "time_infected" <- as.double(NA)
    g %v% "time_recovered" <- as.double(NA)
    g %v% "time_mutated" <- as.double(NA)
    g %v% "infection_history" <- as.integer(NA)
    g %v% "infector_id" <- as.integer(NA)
    g %v% "status" <- "susceptible"
    g %v% "phylo_id" <- as.integer(NA)
    result.nets[[i]] <- .Call("within_epiSimSSA_R", g, trans, infectious.period, 
                              mutation.rate, max.steps, as.integer(verbose), as.double(immune.decay), 
                              as.integer(make.time.series), as.integer(initial.infectious.id), as.double(eff),
                              as.double(conv), as.integer(inoc.size))
  }
  if (i > 1) 
    result.nets
  else result.nets[[1]]
}
