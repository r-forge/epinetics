`net.simpsons.index` <-
function(x) { tot<-sum(table(as.numeric(na.omit(x %v% "infection_history"))));  sum((table(as.numeric(na.omit(x %v% "infection_history")))/tot )^2)}

