full.sc.epi.size <-
function (net, ...) 
{
    p <- strain.degree.dist(net)
    q <- strain.degree.pmf(p)
    u.T.list <- sc.excess.degree.fixed.point(q, ...)
    scd <- sc.dist(net)
    frac <- final.sc.epi.size(u.T.list, q, scd)
    size <- network.size(net)
    frac * size
}
