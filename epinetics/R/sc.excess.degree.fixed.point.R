sc.excess.degree.fixed.point <-
function (tab.list, tol = 0.001, max.iter = 100, base.T = 0.3, 
    inf.decay = 0.1, sus.decay = 0, uninfected.dist = 115, reference.sc = 1) 
{
    num.codes <- length(tab.list)
    codes <- 1:num.codes
    codes[num.codes] <- -uninfected.dist
    u <- array(0.01, dim = rep(num.codes, 2))
    T <- array(0, dim = rep(num.codes, 2))
    for (i in 1:num.codes) {
        for (j in 1:num.codes) {
            infector.dist <- abs(reference.sc - codes[j])
            target.dist <- abs(reference.sc - codes[i])
            T[i, j] <- base.T * (1 - exp(-sus.decay * infector.dist - 
                inf.decay * target.dist))
        }
    }
    print(T)
    u.new <- u
    iter <- 0
    diff <- 1
    while (tol < diff && iter <= max.iter) {
        iter <- iter + 1
        for (i in 1:num.codes) {
            for (j in 1:num.codes) {
                u.new[i, j] <- excess.sc.degree.pgf(u[, i], T[, 
                  i], tab.list[[i]], j)
            }
        }
        diff <- sum(abs(u.new - u))
        print(diff)
        print(u)
        u <- u.new
    }
    list(u, T)
}
