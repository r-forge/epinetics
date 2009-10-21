`get.times.mutated` <-
function(x) rapply(lapply(x$val, function(x) length(x$"time_mutated")-1), function(x) x)

