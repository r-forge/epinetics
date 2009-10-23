.packageName <- 'epineticseq'

.First.lib <-
 function (lib.loc, pkg){
    library.dynam('epineticseq', pkg, lib.loc)
    
    require(sna);
    require(degreenet);
    require(network);
}

.Last.lib <-
 function (lib.loc, pkg){
  library.dynam.unload('epineticseq', pkg, lib.loc)
    
  pos <- match ("package::sna", search())
  if (!is.na(pos)){
    detach (pos = pos)
  }
  pos <- match ("package::degreenet", search())
  if (!is.na(pos)){
    detach (pos = pos)
  }
  pos <- match ("package::network", search())
  if (!is.na(pos)){
    detach (pos = pos)
  }
    
}
