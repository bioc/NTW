`P.preestimation` <-
function(X,topK)
{

    nExp   <- ncol(X)
    nGenes <- nrow(X)
    IX   <- matrix(0,nrow=nGenes,ncol=nExp)
    for(j in 1:nExp)
   {
      colX <- abs(X[,j])
	    srt  <- sort(colX,decreasing=TRUE,index.return=TRUE)
	    IX[srt$ix[1:topK],j] <- 1
   }
   IX
}

