within.param.loop <-
  function (net, reps = 1, ays = 0.00, nus = 1/7, transs = 2/14, 
            thetas = 0.0, mus = 0.001, verbose = 1, conv=.1,
            eff=.1, inoc.size=100, df.out = 1, ...) 
{
  mean.degree <- mean(degree(net, cmode = "indegree"))
  results <- list()
  y <- 0
  counter <- 0
  num.sims <- reps * length(ays) * length(nus) * length(transs) * 
    length(thetas) * length(mus) * length(conv) * length(eff) * length(inoc.size)
  A <- data.frame(a = double(num.sims), nu = double(num.sims), 
                  trans = double(num.sims), theta = double(num.sims), mu = double(num.sims), 
                  es = double(num.sims), length = double(num.sims), max.D = double(num.sims), 
                  mean.D = double(num.sims), var.D = double(num.sims), 
                  simp.recip = double(num.sims), conv=double(num.sims), eff=double(num.sims),
                  inoc=double(num.sims), time.half.peak=double(num.sims), time.peak=double(num.sims),
                  time.half.post.peak=double(num.sims), time.inoc.post.peak=double(num.sims),
                  num.inf.half.peak=integer(num.sims), num.inf.peak=integer(num.sims),
                  num.inf.half.post.peak=integer(num.sims), num.inf.inoc.post.peak=integer(num.sims),
                  frac.ecp.half.peak=double(num.sims), frac.ecp.peak=double(num.sims),
                  frac.ecp.post.peak=double(num.sims), frac.ecp.inoc.post.peak=double(num.sims),
                  total.ecp=integer(num.sims))

  
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
#            results[[counter]] <- within.epi.sim.ssa(net, reps = reps, 
#                                              trans = trans, max.steps = 1e+05, infectious.period = 1/nu, 
#                                              mutation.rate = mu, verbose = 0, immune.decay = a, 
#                                              initial.susceptibility = (1 - theta), inoc.size=inoc.size,
#                                                     eff=eff, conv=conv)
            results[[counter]] <- within.epi.sim.ssa(net=net, trans=trans, mu=mu, eff=eff, conv=conv,
                                                     reps=reps, infectious.period=1/nu, inoc.size=inoc.size,... )
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
                  inoc <- sum(na.omit(g %v% "time_infected" == 0))
                  mts <- make.multi.phy.time.series(g)
                  es <- sum(g %v% "status" %in% "recovered")
                  A[row.num, ] <- c(a, nu, trans, theta,
                                    mu, es, 
                                    max(na.omit(g %v% "time_recovered")), 
                                    max(na.omit(g %v% "infection_history")), 
                                    mean(na.omit(g %v% "infection_history")), 
                                    var(na.omit(g %v% "infection_history")), 
                                    1/net.simpsons.index(g),
                                    g %n% "conversion_factor",
                                    g %n% "efficiency",
                                    inoc,
                                    mts.cruncher(mts, inoc),
                                    sum(na.omit(g %v% "phylo_id" != 1))/es
                                    )

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
