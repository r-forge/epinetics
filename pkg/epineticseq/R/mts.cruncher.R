mts.cruncher <- function(mts, inoc){
 
 max.all <- max(mts[,2])
 t.peak <- double(1)
 num.peak <- integer(1)
 frac.peak <- as.double(0)

 t.half.peak <- double(1)
 num.half.peak <- integer(1)
 frac.half.peak <- as.double(0)

 t.half.post.peak <- double(1)
 num.half.post.peak <- integer(1)
 frac.half.post.peak <- as.double(0)

 t.inoc.post.peak <- double(1)
 num.inoc.post.peak <- integer(1)
 frac.inoc.post.peak <- as.double(0)
  
 flag1 <- 1
 flag2 <- 0
 
 apply(mts, 1, function (x)
        {
         if(x[2]==max.all)
           {
             t.peak<<-x[1];
             num.peak<<-x[2];
             frac.peak<<-(x[2]-x[3])/x[2];
           }
         else if (x[2] > max.all/2 && flag1 == 1)
           {
             t.half.peak<<-x[1];
             num.half.peak <<-x[2];
             frac.half.peak<<-(x[2]-x[3])/x[2];
             flag1 <<- 0
             flag2 <<- 1
           }
         else if (x[2] < max.all/2 && flag2 == 1)
           {
             t.half.post.peak<<-x[1];
             num.half.post.peak<<-x[2];
             frac.half.post.peak<<-(x[2]-x[3])/x[2]
             flag2 <<-0
           }
         else if (x[2] <= inoc && flag1 == 0 && flag2 == 0)
           {
             t.inoc.post.peak<<-x[1];
             num.inoc.post.peak<<-x[2];
             frac.inoc.post.peak<<-(x[2]-x[3])/x[2]
             flag1 <<- 1
           }
        
        })       
 c(t.half.peak, t.peak, t.half.post.peak, t.inoc.post.peak,
   num.half.peak, num.peak, num.half.post.peak, num.inoc.post.peak,
   frac.half.peak, frac.peak, frac.half.post.peak, frac.inoc.post.peak)
}
