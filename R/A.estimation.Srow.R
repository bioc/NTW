`A.estimation.Srow` <-
function(r,cMM.corrected,pred.net,X,P.known,topD,restK,cFlag,sup.drop,noiseLevel)
{
	if(cMM.corrected==1)
	{
		TA<-which(pred.net[r,]>0)
		TA<-sort(union(TA,r))
		}
	else
	{
		TA=1:nrow(X)
		}
	best.rA <- matrix(0,1,nrow(X))
	fun=function(pert,r,X,topD,restk,cFlag,TA,sup.drop,noiseLevel)
	{
		if(sup.drop==-1)
		     return(backward(r,X,pert,topD,restk,cFlag,TA,noiseLevel))
		if(sup.drop==1)
		     return(forward(r,X,pert,topD,restk,cFlag,TA,noiseLevel))
		}
	#restK[r] the restK  of the rth line.
	result <- fun(P.known[r,],r,X,topD,restK[r],cFlag,TA,sup.drop,noiseLevel)
	best.rA <- result$A.row
	return(list(A.row = best.rA))
	}

