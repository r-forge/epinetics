
net.shuffle <-
  function (net, steps=0, metric.interval=0, test.depth.increment=0,
            q.plus = 0.25, test.depth=2, print.metrics=0){
    on.exit(.C("clean_up", PACKAGE="epinetics"));
    edges <- as.matrix(net, matrix.type="edgelist")
    heads <- as.integer(edges[,1])
    tails <- as.integer(edges[,2])
    degree.sum <- as.integer(length(edges))
    size <- as.integer(network.size(net))
    if (steps == 0 || steps < 0)
      num.steps = as.integer (degree.sum * 50)
    else num.steps = as.integer(steps)
    if (metric.interval == 0)
      metric.interval = as.integer(num.steps)
    else metric.interval = as.integer(metric.interval)
   
   swaps.done = 0;
   failure.rate = 1;
   mean.cc = 0;
   assortativity = 0;
   assort.sigma = 0;
   diameter = 0;
   cpl = 1;

    x <- .C("shuffle_net_R", size, degree.sum, heads, tails, 
            num.steps, metric.interval, as.integer(test.depth.increment),
            as.double(q.plus), as.integer(test.depth), as.integer(print.metrics),
	    as.integer(swaps.done), as.double(failure.rate), as.double(mean.cc), as.double(assortativity), 
	    as.double(assort.sigma), as.integer(diameter), as.double(cpl),PACKAGE="epinetics")
    ret <-list()
    g<-network(cbind(x[[3]],x[[4]]), directed=FALSE)
    if (is.connected(g))
      ret[["network"]]<-network(cbind(x[[3]],x[[4]]), directed=FALSE)
    else print ("Error: graph is not connected at end of shuffling\n")   
    ret[["swaps.done"]]<-x[[11]]
    ret[["failure.rate"]]<-x[[12]]
   
if (print.metrics){    
     ret[["mean.cc"]]<-x[[13]]
      ret[["assortativity"]]<-x[[14]]
      ret[["assort.sigma"]]<-x[[15]]
      ret[["diameter"]]<-x[[16]]
      ret[["cpl"]]<-x[[17]]
  }

    ret

}

