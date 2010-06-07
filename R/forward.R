`forward` <-
function(r,X,pert,topD,restk,cFlag,TA,noiseLevel)
{
	no_experiments <- ncol(X)
	no_genes <- nrow(X)
	no_existed=length(TA)

	restk=max(restk-no_existed,0)
	selected_genes <- matrix(TA,1,)
	minCrtValue <- Inf
	bestComb    <- NULL
	if(restk==0)
	{
    if(cFlag=="geo") 
    temp.result<-method.geo(index.vars=TA,X=X,pert=pert)
    if(cFlag=="sse") 
    temp.result<-method.sse(index.vars=TA,X=X,pert=pert)
    if(cFlag=="ml") 
    temp.result<-method.ml(index.vars=TA,X=X,pert=pert,noiseLevel=noiseLevel)
    sol  <- temp.result$sol
    error<- temp.result$error
   	rA=rep(0,no_genes)
   	rA[TA] <- sol
   	return(list(A.row=rA, CrtValue=error))
    } else
    {
   	  for (k in 1:restk)
    	{
    		comb1 <- NULL
    		for (g in 1:nrow(selected_genes))
    		{
    			rest_genes <- setdiff(1:no_genes,selected_genes[g,])
    			if (length(rest_genes)==1)
    			{
			      comb <- rest_genes
          } else
    			{
			      comb  <- t( combn(rest_genes,1))
          }
    			comb_temp <- comb.matrix(selected_genes[g,],comb)
    			comb1 <- rbind(comb1,comb_temp)
    			}
    		comb = unique(t(apply(comb1,1,sort)))
    		sol <- matrix(0,nrow = nrow(comb),ncol=ncol(comb))
    		error <- matrix(0,nrow = nrow(comb),1)
    		for (n in 1:nrow(comb))
    		{
    
          if(cFlag=="geo") 
          temp.result<-method.geo(index.vars=comb[n,],X=X,pert=pert)
          if(cFlag=="sse") 
          temp.result<-method.sse(index.vars=comb[n,],X=X,pert=pert)
          if(cFlag=="ml") 
          temp.result<-method.ml(index.vars=comb[n,],X=X,pert=pert,noiseLevel=noiseLevel)
          sol[n,]<-  temp.result$sol
          error[n,]<-temp.result$error
        }
    		srtE <- sort(error,index.return=TRUE)
    		topd=topD
    		topD=min(topD,nrow(comb))
    		topD_comb <- comb[srtE$ix[1:topD],]
    		topD_sol=sol[srtE$ix[1:topD],]
    		topD_error=error[srtE$ix[1:topD],]
    		if (is.vector(topD_comb))
    		{
    			topD_comb <- t(as.matrix(topD_comb))
    			topD_sol <- t(as.matrix(topD_sol))
   			}
    		topD_error <- t(as.matrix(topD_error))
    		selected_genes <-topD_comb
    		crtValue <- srtE$x[1]
    		# select the combination best fit for the criterion
    		if(minCrtValue > crtValue)
    		{
    			minCrtValue <- crtValue
    			bestComb    <- list(selected_genes[1,],topD_sol[1,],topD_error[1,])
   			}
    		topD=topd
      }
    	rA=rep(0,no_genes)
	   	if(minCrtValue!=Inf)
	    {
      	rA[as.vector(bestComb[[1]])]=as.vector(bestComb[[2]])
      }
    	return(list(A.row=rA, CrtValue=minCrtValue))
    }
}

