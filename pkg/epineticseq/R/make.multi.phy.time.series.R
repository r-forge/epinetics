make.multi.phy.time.series <-
  function (sim.result, time.steps = 100) 
{
  all.ts <- make.phy.time.series(sim.result)
  wt.ts <- make.phy.time.series(sim.result, all.phy=FALSE, phylo.vector=1)
  res <- cbind(all.ts, wt.ts[,-1])
  colnames(res) <- c("t", "all.inf", "wt.inf")
  res
}
