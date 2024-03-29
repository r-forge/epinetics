
cond.hist <-
  function(A){

require(lattice)

sundar.theme <- function() {
        par <- col.whitebg()
        par$strip.background$col <- rep("#000099", 7)
        par$add.text$col <- "#eeeeaa"
        par$add.text$font <- 2
        par$background$col <- "#ffffff"
        par$superpose.line$lty <- rep(1, 7)
        par$superpose.line$col[1:2] <- c("#880000", "#008800")
        par$superpose.symbol$col[1:2] <- c("#880000", "#008800")
        par
      }

trellis.par.set(sundar.theme())

print(  # necessary if the file is source()'d
  histogram( ~ mean.D | paste("mu=",as.character(mu),sep="") + paste("T=",as.character(round(trans,2)),sep=""), data = A,
          xlab ="mean drift", type = "percent",
          panel = function(x, ...) {
              panel.histogram(x, ...)
#              panel.mathdensity(dmath = dnorm, col = "black",
#                                args = list(mean=mean(x),sd=sd(x)))
          } )
  )

}
