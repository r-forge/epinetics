`vt.degree.pgf` <-
  function (network, x, T, nu=.2){
    deg.set <- degree(network, cmode="indegree")
    deg.dist <- table(deg.set)
    deg.dist.norm <- deg.dist/sum(deg.dist)
    beta <- beta.from.T(T, nu)
    integrand <- function(T.real){
      inf.period.real <- log(1-T.real)/-beta
      prob.T.real = nu*exp(-nu*inf.period.real)
      sum <- 0
      for(i in 1:length(deg.dist.norm)){
        k <- as.numeric(names(deg.dist.norm)[i])
        pk <- as.numeric(deg.dist.norm[i])
        sum <- sum + prob.T.real*pk*(1+(x-1)*T.real)^(k)
      }
      sum
    }
    integrate(integrand, lower=0.001, upper=1)
  }

