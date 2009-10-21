param.loop <-
  function (net, reps = 1, ays = 0.01, nus = 0.2, transs = 0.24, 
            thetas = 0.1, mus = 0.05, verbose = 1, df.out = 1) 
{
  mean.degree <- mean(degree(net, cmode = "indegree"))
  results <- list()
  y <- 0
  counter <- 0
  num.sims <- reps * length(ays) * length(nus) * length(transs) * 
    length(thetas) * length(mus)
  A <- data.frame(a = double(num.sims), nu = double(num.sims), 
                  trans = double(num.sims), theta = double(num.sims), mu = double(num.sims), 
                  es = double(num.sims), length = double(num.sims), max.D = double(num.sims), 
                  mean.D = double(num.sims), var.D = double(num.sims), 
                  simp.recip = double(num.sims))
  for (i in 1:length(ays)) {
    a <- ays[i]
    for (j in 1:length(nus)) {
      nu <- nus[j]
      for (k in 1:length(transs)) {
        trans <- transs[k]
        if (trans >= 1) 
          trans <- 0.99
        for (l in 1:length(thetas)) {
          theta <- thetas[l]
          for (m in 1:length(mus)) {
            mu <- mus[m]
            if (verbose) 
              print(paste("a:", a, " nu:", nu, " trans:~", 
                          round(trans, 2), " theta:", theta, " mu:", 
                          mu))
            counter = counter + 1
            results[[counter]] <- epi.sim.ssa(net, reps = reps, 
                                              trans = trans, max.steps = 1e+05, infectious.period = 1/nu, 
                                              mutation.rate = mu, verbose = 0, immune.decay = a, 
                                              initial.susceptibility = (1 - theta))
            if (df.out) {
              for (r in 1:counter) {
                for (p in 1:reps) {
                  if (reps > 1) {
                    g <- results[[r]][[p]]
                  }
                  else {
                    g <- results[[r]]
                  }
                  row.num <- y + r + (p - 1)
                  A[row.num, ] <- c(a, nu, trans, theta, 
                                    mu, sum(g %v% "status" == "recovered"), 
                                    max(na.omit(g %v% "time_recovered")), 
                                    max(na.omit(g %v% "infection_history")), 
                                    mean(na.omit(g %v% "infection_history")), 
                                    var(na.omit(g %v% "infection_history")), 
                                    1/net.simpsons.index(g))
                }
              }
              y <- y + (p - 1) + r
              counter = 0
            }
          }
        }
      }
    }
  }
  print(paste("num sims:", num.sims))
  if (df.out) {
    A
  }
  else results
}
