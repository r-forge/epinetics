strain.degree.pgf <-
  function (z, T, tab) 
{
  if (is.null(tab)) 
    stop("Pmf table is null. Does any node have this strain code?\n")
  num.codes <- length(dimnames(tab))
  if (num.codes != length(z)) 
    stop("Length of vector argument z must equal the dimensions of the pgf\n")
  z.exp <- vector(mode = "numeric", length = num.codes)
  dummy <- 1 + z * T - T
  ft <- ftable(tab)
  col.vals <- as.numeric(attributes(ft)$col.vars[[1]])
  index <- 1
  sum <- 0
  row.mat <- expand.grid(lapply(attributes(ft)$row.vars, as.numeric))
  if (dim(row.mat)[2] > 1) {
    row.mat <- row.mat[do.call(order, row.mat), ]
  }
  rows.in.row.mat <- dim(row.mat)[1]
  for (i in 1:length(col.vals)) {
    z.exp[num.codes] <- col.vals[i]
    for (j in 1:rows.in.row.mat) {
      z.exp[-num.codes] <- as.numeric(row.mat[j, ])
      sum <- sum + ft[index] * prod(dummy^z.exp)
      index <- index + 1
    }
  }
  sum
}
