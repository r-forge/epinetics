dist.at.time <-
  function(net, time=9.7) {
    for (i in 1:length(net %n% "time_series"))
      {
        t <- (net %n% "time_series")[i]
        if(t >= time) break; t <- NA
      }
    if(!(is.na(t)))
      {
        dist <- (net %n% "mean_distance_series")[i]
        dist
      }
    else NA
  }

