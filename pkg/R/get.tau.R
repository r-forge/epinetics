get.tau <-
  function(trans=0.5, nu=0.2, kappa=14*.96, N=1000){
    beta <- beta.from.T(T=trans, nu=nu)
    alpha <- beta*kappa - nu
    log(sqrt(N))/alpha
  }
