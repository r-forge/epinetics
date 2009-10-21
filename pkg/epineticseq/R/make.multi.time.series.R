make.multi.time.series <-
function (sim.result, time.steps = 100) 
{
    max.strain.id <- max(na.omit(sim.result %v% "infection_history"))
    min.strain.id <- min(na.omit(sim.result %v% "infection_history"))
    mts <- make.time.series(sim.result, all.strains = FALSE, 
        strain.vector = min.strain.id, time.steps = time.steps)
    vec <- (min.strain.id + 1):max.strain.id
    for (id in vec) {
        next.ts <- make.time.series(sim.result, all.strains = FALSE, 
            strain.vector = id, time.steps = time.steps)
        mts <- cbind(mts, next.ts[, 2])
    }
    colnames(mts) <- c("time", paste("I", min.strain.id:max.strain.id, 
        sep = "."))
    mts
}
