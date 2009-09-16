mean.sc.degree <-
  function (tab, sc) 
{
  if(is.null(tab))
    stop("Pmf table is null. Does any node have this strain code?\n")
  mtab <- margin.table(tab, sc)
  sum <- 0
  for (i in 1:length(mtab)) {
    sum <- sum + as.numeric(names(mtab)[i]) * mtab[i]
  }
  names(sum) <- NULL
  sum
}
