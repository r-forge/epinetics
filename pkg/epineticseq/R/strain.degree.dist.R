`strain.degree.dist` <-
function (net) 
{
  N <- network.size(net)
  types <- unique(net %v% "infection_history")
  default.strain <- (net %n% "this_years_strain_0_was_last_years")
  max.types <- max(na.omit(types))
  for (i in 1:N) {
    hood <- get.neighborhood(net, i)
    hood.strains <- (net %v% "infection_history")[hood]
    hood.strains <- hood.strains + 1
    num.uninfected <- length(hood.strains[hood.strains %in% 
                                          NA])
    degrees <- tabulate(na.omit(hood.strains), nbins = max.types + 
                        2)
    degrees[max.types + 2] <- num.uninfected
#    print(degrees)
    net$val[[i]]$strain_degrees <- degrees
  }
  print("Degrees types start at strain 0 and count up to the highest strain\n and then the number of uninfected neighbors is appended on the end.")
  net
}

