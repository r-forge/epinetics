custplot <- function(a){
 pdf()
# plot(a$t, a$"cases", col=4, lty=1)
 plot(a$t, a$D, xlab="time (days)", ylab="", type="l")
 lines(a$t, a$fs, col=3)
 lines(a$t, a$e.D, col=2)
 dev.off()
}
