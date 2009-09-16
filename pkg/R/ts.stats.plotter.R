`ts.stats.plotter` <-
function (mts, mu=0.05){
    stats2 <- mts.stats(mts)
    plot(stats2[,1],stats2[,1]*mu, xlab="days", ylab="")
    points(stats2[-c(1,101),1],stats2[-c(1,101),2], pch="v")
    points(stats2[-c(1,101),1],stats2[-c(1,101),3], pch="m")
  }

