`AP.estimation.Srow` <-
function(r,cMM.corrected,pred.net,X,IX,topD,restK,cFlag,sup.drop,numP,noiseLevel)
{
	#print(numP)
	if(cMM.corrected==1)
  {
		#vec: row vector
		TA<-which(pred.net[r,]>0)
		##Change TA to include r
        TA<-sort(union(TA,r))###why r should be included?
  } 
  else
  {
    TA=1:nrow(X)
  }

	nExps   <- ncol(X)
	best.rA <- matrix(0,1,nrow(X))
	minCrtValue <- Inf
	idx4cP <- NULL

	fun=function(pert,r,X,topD,restk,cFlag,TA,sup.drop,noiseLevel)
	{
		if(sup.drop==-1)
		     return(backward(r,X,pert,topD,restk,cFlag,TA,noiseLevel))
		if(sup.drop==1)
		     return(forward(r,X,pert,topD,restk,cFlag,TA,noiseLevel))
	}

	#If one row of IX is all zero, then the corresponding row of P matrix is all zero
	vec1=IX[r,]
	
	if(length(vec1[vec1>0])==0)
	   idx4cP<-NULL
	#If one row of IX is not all zero
	else
	{
		#vec4cP: the set of all the none-zero places of the rth row of P matrix
		vec4cP<-which(vec1>0)
		len<-length(vec4cP)
		if(numP>len)
		    numP<-len
		for(j in 1:numP)
		{
			if(len!=1)
			{
				matrixP<-t(combn(vec4cP,j))
				}
			else
			{
				matrixP<-as.matrix(vec4cP)
				}
			#restK[r] the restK  of the rth line.
			# estimate a single row of A with P fixed
			#Compute matrixPfull
			matrixPfull<-matrix(0,nrow=nrow(matrixP),ncol=ncol(X))
			for(k in 1:nrow(matrixPfull))
			{
				matrixPfull[k,matrixP[k]]=1
				}
			result=apply(matrixPfull,MARGIN=1,fun,r,X,topD,restK[r],cFlag,TA,sup.drop,noiseLevel)
			for(i in 1:length(result))
			{
				if(result[[i]]$CrtValue < minCrtValue)
				{
					best.rA <- result[[i]]$A.row
					minCrtValue <- result[[i]]$CrtValue
					idx4cP <- matrixP[i,]
					}
				}
			}
		}
	return(list(P.index=idx4cP,A.row = best.rA))
	}

