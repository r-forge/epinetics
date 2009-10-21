seq.list.to.nexus <-
  function(isolation.data, sequence.list, filename="simulated.nex"){
    cat("#NEXUS\n\n", file=filename)
    cat("BEGIN DATA;\n", file=filename, append=TRUE)
    tab <- "        "
    cat(tab, "DIMENSIONS NTAX=", length(isolation.data),
        " NCHAR=", length(sequence.list[[1]]), ";\n", file=filename, sep="", append=TRUE)
    cat(tab, "FORMAT MISSING=? GAP=- DATATYPE=DNA;\n", file=filename, sep="", append=TRUE)
    cat(tab, "MATRIX\n", file=filename, sep="", append=TRUE)
    for (i in 1:length(isolation.data)){
      cat(tab, "host", isolation.data[[i]][1], "_time", isolation.data[[i]][2], tab, sequence.list[[isolation.data[[i]][3]]], "\n", file=filename, sep="", append=TRUE)
    }
    cat(";\nEnd;\n", file=filename, sep="", append=TRUE)
  }
    
    
    
