`NTW` <-
function(X,restK,topD,topK=NULL,P.known=NULL,cFlag,pred.net=NULL,sup.drop=-1,numP=NULL,noiseLevel=0.1)
{
	
  if (!is.matrix(X))
       stop("X must be a matrix!\n")
       
  if(length(restK)!=nrow(X))
    stop("The length of restK should equal to the number of rows of the X matrix.")
    
       
  if((cFlag!="geo")&&(cFlag!="sse")&&(cFlag!="ml"))
    stop("An invalid cFlag indicator. cFlag should be one of 'geo', 'sse' and 'ml'.")
  
  if(is.null(pred.net))
  {
  	if(sup.drop==1)
  	   stop("Since prior gene association information is not available (pred.net=NULL), a backward method should be used as default (sup.drop = -1)")
  	}
  else
  {
  	if(nrow(pred.net)!=nrow(X) || ncol(pred.net)!=nrow(X))
  	    stop("pred.net is the prior gene-gene interaction network. it's a square matrix, of which the dimension equals to the row length of X. ")
  	}
  	
   if(!is.null(P.known))
  {
  	if(nrow(P.known)!= nrow(X) || ncol(P.known)!= ncol(X))
      stop("P.known should has the same dimensions as X.")
  	
      if(!is.null(topK)||!is.null(numP))
     warning("topK or numP is useless when P matrix is known.")
  	}
  else
  {
  	if(is.null(topK)||is.null(numP))
  	  stop("Invalid input for topK or numP. When P matrix is required for estimation, topK and numP should be set as an integer.")
  	
  	} 
  
  
  if(max(restK)>nrow(X))
  {
  	warning("At least one of the element in restK is larger than the gene counts (nrow(X))! The element, which is larger than nrow(X), in restK is useless for regression.")
  	
  	}
  if(!is.null(topK))
  if(topK>nrow(X))
  {
  	warning("topK should be smaller than the gene counts (nrow(X)).")
  	
  	}
      
  if(is.null(pred.net))
    cMM.corrected=0 
  else
    cMM.corrected=1
  if(is.null(P.known))
  {
    est.A=matrix(0,nrow(X),nrow(X))
    est.P=matrix(0,nrow(X),ncol(X))
    IX <- P.preestimation(X,topK)
    #Compute all rows in X
    mrow=as.matrix(1:nrow(X))
    #MARGIN=1 by rows
    est=apply(mrow,MARGIN=1,FUN=AP.estimation.Srow,cMM.corrected=cMM.corrected,pred.net=pred.net,X,IX=IX,topD=topD,restK=restK,cFlag=cFlag,sup.drop=sup.drop,numP=numP,noiseLevel=noiseLevel)
    for(i in 1:nrow(X))
    {
      est.A[i,]=est[[i]]$A.row
      if(length(est[[i]]$P.index)!=0)
        est.P[i,est[[i]]$P.index]=1
    }
    #set the row name and column name of A matrix and P matrix
    dimnames(est.A)=list(rownames(X),rownames(X))
    dimnames(est.P)=dimnames(X)
    return(list(est.A=est.A,est.P=est.P))
  } 
  else
  {
      est.A=matrix(0,nrow(X),nrow(X))
      #Compute all rows in X
      mrow=as.matrix(1:nrow(X))
      #MARGIN=1 by rows
      est=apply(mrow,MARGIN=1,FUN=A.estimation.Srow,cMM.corrected=cMM.corrected,pred.net=pred.net,X,P.known=P.known,topD=topD,restK=restK,cFlag=cFlag,sup.drop=sup.drop,noiseLevel=noiseLevel)
      for(i in 1:nrow(X))
      {
        est.A[i,]=est[[i]]$A.row
      }
      #set the row name and column name of A matrix 
      dimnames(est.A)=list(rownames(X),rownames(X))
      return(list(est.A=est.A))    	    
    }
  
}

