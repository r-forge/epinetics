final.epi.size <-
function (network, T, tol = 0.001, max.iterations = 10) 
{
    u <- 0
    diff <- 1
    counter <- 0
    while (diff > tol && counter < max.iterations) {
        u.new <- excess.degree.pgf(network = network, T = T, 
            x = u)
        diff <- abs(u.new - u)
        u <- u.new
        counter <- counter + 1
    }
    print(u)
    network.size(network) * (1 - degree.pgf(network = network, 
        x = u, T = T))
}
