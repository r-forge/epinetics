simulate.sequences <-
  function(phylo.tree,
           ancestral.sequence,
           rate.mat=matrix(c(
             c(0,1,1,1),
             c(1,0,1,1),
             c(1,1,0,1),
             c(1,1,1,0)),4,4,
             dimnames=list(
               c("a","g","t","c"),
               c("a","g","t","c")))){
    num.nodes <- network.size(phylo.tree)
    root.id <- which(is.na(as.numeric(phylo.tree$iel)))
    sequence.list <<- list()
    sequence.list[[root.id]] <<- ancestral.sequence
    generated <- vector(mode="logical", length=num.nodes)
    generated[root.id] <- TRUE
    flag <- 0
    sequence.bot(root.id, phylo.tree, generated, rate.mat, flag)
    sequence.list
  }
    
    
  

