sequence.bot <-
  function(node.id, phylo.network, generated, rate.mat, flag){
    parent.id <- get.neighborhood(x=phylo.network, v=node.id, type="in")
    if (!generated[node.id]){
      sequence.list[[node.id]] <<- sequence.list[[parent.id]]
      mutate.sequence(node.id, rate.mat)
      generated[node.id] <- TRUE 
    }
    child.nodes <- get.neighborhood(x=phylo.network,v=node.id, type="out")
    if(length(child.nodes)){
      for (i in 1:length(child.nodes)){
        x <- child.nodes[i]
        if (!generated[x]){
          sequence.bot(x, phylo.network, generated, rate.mat, flag)
          return(0)
        }
      }

    }
    if (length(parent.id)!=0){
      sequence.bot(parent.id, phylo.network, generated=generated, rate.mat, flag)
    }
  }
      
      
