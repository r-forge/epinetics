beta.from.T <- function(T=.2, nu=0.2){
  T*nu/(1-T)
}
E.v <-
  function(b.c=2.8, nu=0.2, t=3){
    a <- b.c-nu
    prod <- 1/(2*b.c/(nu+b.c)-1)
    prod <- prod*exp(a*t)
    prod <- prod*(exp(a*t)-1)
    prod
  }
`i.of.t` <-
  function (sim.result, times = 10, all.strains = TRUE, strain.vector = NA) 
{
  time.sequence <- times
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

make.inf.tau.vec <-
  function(net.list, tau){
     rapply(lapply(net.list, function(x) i.of.t(x ,times=tau)[1,2]), function(x) x)
   }
