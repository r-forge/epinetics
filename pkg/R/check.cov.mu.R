check.cov.mu <-
  function(net){
    mu <- net %n% "mutation_rate"
    df <- data.frame(t=net%n%"time_series", cov=net%n%"cov_series", D=net%n%"mean_distance_series", fs=net%n%"fitness_series", cases=net%n%"case_series")
    len <- dim(df)[1]
    df$del.t <- -1*c(df$t[-len]-df$t[-1],NA)
    df$del.D <- -1*c(df$D[-len]-df$D[-1],NA)
    df$mu.del.t <- mu * df$del.t
    df$cov.del.t <- df$cov * df$del.t
    df$e.del.D.del.t <- df$mu.del.t + df$cov.del.t
    df$error <- df$e.del.D.del.t - df$del.D
    df$e.D <- df$e.del.D.del.t
    for (i in 2:length(df$e.D)){
      df$e.D[i] <- df$e.del.D.del.t[i]+df$e.D[i-1]
    }
    print(mean(na.omit(df$error)))
    df
  }
