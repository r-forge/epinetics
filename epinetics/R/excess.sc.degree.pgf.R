excess.sc.degree.pgf <-
  function(z, T, tab, sc){
  if(is.null(tab)){
#    stop("Pmf table is null. Does any node have this strain code?\n")
# return 1 in this case since if they are not in the network they can be considered as
    #having excess degree 0 with probability 1, and anything^0 = 1 
    return(1)
  }
  num.codes <- length(dimnames(tab))
  if(num.codes!=length(z)){
    stop("Length of vector argument z must equal the dimensions of the pgf\n")
  }
  z.exp <- vector(mode = "numeric", length = num.codes)
  if(sc > num.codes || sc < 1)
    stop("Strain code (argument sc) must be in [1,num.codes]\n")

  denom <- mean.sc.degree(tab, sc)
  if (denom==0) return(1)
  dummy <- 1 +z*T - T
  ft <- ftable(tab)
#  print(ft)
  col.vals <- as.numeric(attributes(ft)$col.vars[[1]])
  index <- 1
  sum <- 0
  row.mat <- expand.grid(lapply(attributes(ft)$row.vars, as.numeric))
  if(dim(row.mat)[2]>1){
    row.mat <- row.mat[do.call(order, row.mat), ]
  }
#  print(row.mat)
  rows.in.row.mat <- dim(row.mat)[1]
  for (i in 1:length(col.vals)) {
    flag <- 0
    if(sc==num.codes){
      weight <- col.vals[i]/denom
      if(denom==0) weight <- 0
      z.exp[num.codes] <- col.vals[i]-1
      flag <- 1
    }
    else{
      z.exp[num.codes] <- col.vals[i]
    }
    for (j in 1:rows.in.row.mat) {
      z.exp[-num.codes] <- as.numeric(row.mat[j, ])
      if(flag==0){
        weight <- z.exp[sc]/denom
        if(denom==0) weight <- 0
        z.exp[sc] <- z.exp[sc] - 1
      }
#      print(z.exp)
#      print(ft[index])
      # only add to sum if excess degrees are at least 0 to avoid infinite terms
      if(all(z.exp >= 0))
        sum <- sum + weight * ft[index] * prod(dummy^z.exp)
      index <- index + 1
#      print(z.exp)
#      print(c(sum, weight, ft[index], prod(z^z.exp)))
    }
  }
  sum
}
