`immunize` <-
function (ego.id, strain.id){
    set.vertex.attribute (epi.network, "immune_memory", strain.id, ego.id)
  }

