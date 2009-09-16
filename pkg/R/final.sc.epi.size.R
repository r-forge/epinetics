final.sc.epi.size <-
function (u.T.list, tab.list, sc.dist) 
{
    u <- u.T.list[[1]]
    T <- u.T.list[[2]]
    num.codes <- length(tab.list)
    sum <- 0
    for (i in 1:num.codes) {
        if (!is.null(tab.list[[i]])) {
            sum <- sum + as.numeric(sc.dist[i]) * strain.degree.pgf(u[, 
                i], T[, i], tab.list[[i]])
        }
    }
    1 - sum
}
