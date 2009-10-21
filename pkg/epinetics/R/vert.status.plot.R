vert.status.plot <-
   function (net){
   plot(net, displaylabels=TRUE, vertex.col=match(net %v% "status", levels(as.factor((net %v% "status"))) ))
}
