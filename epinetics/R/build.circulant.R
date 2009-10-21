build.circulant <-
  function (N=10, mean.k=2){
    m <- matrix(0, N, N)
    K <- mean.k / 2
    if (mean.k %% 2 )    #if odd
      print(paste("warning: only even mean degrees supported"))
    for (j in 1:N){
      nb = K
      non.nb = N - 1 - 2*K
      if (j + 1 <= N){
        for ( i in (j+1):N){
          if (nb){
            m[i,j]=1
            nb = nb -1
          }
          else if (non.nb){
            m[i,j]=0
            non.nb = non.nb - 1
          }
          else{
            m[i,j]=1
          }
        }
      }
      if (j-1 >= 1)
        for ( i in 1:(j-1)){
          if (nb){
            m[i,j]=1
            nb = nb -1
          }
          else if (non.nb){
            m[i,j]=0
            non.nb = non.nb - 1
          }
          else{
            m[i,j]=1
          }
      }
   }
   if (nb || non.nb) print("error in indexing")
    else network(m, directed=FALSE)
  }
                 

