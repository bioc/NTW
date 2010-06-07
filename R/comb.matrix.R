`comb.matrix` <-
function(x,y)
{
# Create all combinations of vectors especially for matrices
# The resultant matrix is nrow(x)*nrow(y) rows and ncol(x)+ncol(y) columns
   if(!is.matrix(x))
      x <- t(as.matrix(x))
   if(!is.matrix(y))
      y <- t(as.matrix(y))
   r1 <- nrow(x)
   r2 <- nrow(y)

   tmp<- NULL
# x varys fastest
   for (i in 1:r2){
       tmp <- rbind(tmp,x)
   }
# each row of y should repeat r1 times
   tmp2 <- NULL
   for (i in 1:r2){
       for (j in 1:r1){
            tmp2 <- rbind(tmp2,y[i,])
       }
   }
  res <- cbind(tmp,tmp2)
  res
}
