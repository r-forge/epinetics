`cust3plot` <-
function (net, marginal.hist = 1, cov.line = 0, ret.df=1) 
{
  def.par <- par(no.readonly = TRUE)
  df <- check.cov.mu(net)
  dist <- as.numeric(na.omit(net %v% "infection_history"))
  par(fig = c(0, 0.8, 0, 0.8))
  if (marginal.hist == 1) 
    ylims <- c(0, max(dist))
  else ylims <- c(0, max(df$D))
  plot(df$t, df$D, xlab = "time (days)", ylab = "distance (amino acids)", 
       type = "l", ylim = ylims)
  if (marginal.hist == 1) {
    abline(a = mean(dist), b = 0, lty = 3)
  }
  e.D.0.intercept <- c(0, df$e.D[-length(df$e.D)])
  lines(df$t, e.D.0.intercept, lty = 5)
  par(fig = c(0, 0.8, 0.625, 1), new = TRUE)
  total <- network.size(net)
  frac.inf <- df$cases/total
  if (cov.line == 1) 
    ylims <- c(-1, 1)
  else ylims <- c(0, 1) 
  plot(df$t, df$cases/total, ylim = ylims, xaxt = "n", type = "l", 
       xlab = "", ylab = "frequency")
  if (cov.line == 1) 
    lines(df$t, df$cov, lty = 2)
  dhist <- hist(dist, breaks = seq(min(dist), max(dist), 1), 
                plot = FALSE)
  par(fig = c(0.65, 1, 0, 0.8), new = TRUE)
  if (marginal.hist == 1) {
    barplot(dhist$counts, yaxt = "n", space = 0, horiz = TRUE, 
            xlab = "total infected")
  }
  par(def.par)
  if(ret.df)
    df
}

