sample.sequences <-
  function(epidemic.model, sample.size){
    infecteds <- which(epidemic.model %v% "status" %in% "recovered")
    if (sample.size > length(infecteds)){
      stop("Sampling more hosts than were infected\n")
    }
    x <- sample(infecteds, sample.size, replace=FALSE)
    data <- list()
    vertex.attribs <- epidemic.model$val
    for( i in 1:length(x)){
      host.id <- x[i]
      isolation.time <- runif(1,min=vertex.attribs[[host.id]]$time_infected,
                              max=vertex.attribs[[host.id]]$time_recovered)
      y <- 2
      while(isolation.time > vertex.attribs[[host.id]]$time_mutated[y] && y <= length(vertex.attribs[[host.id]]$time_mutated)){
        y <- y + 1
      }
      data[[i]] <- c(host.id, isolation.time, vertex.attribs[[host.id]]$phylo_id[y-1])
    }
    data
}
