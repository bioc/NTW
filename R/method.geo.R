`method.geo` <-
function(index.vars,X,pert)
{
	opt=function(bta)
	{
		if(length(bta[bta==0])!=0)
		{
			return(100000)
			break
		} else
		{
      #geometrix mean
		  res <-as.vector(-pert) - as.vector(matrix(bta, 1, ) %*% r1)
			geomean<-sum(res^2)/(abs(prod(bta))^(2/length(bta)))
			return(geomean)
	  }
  }
	r1=X[index.vars,]
	rcd=1/kappa(r1%*%t(r1))
	if(rcd>=10^(-9))
	{
	  start_value=-pert %*% t(r1) %*% solve(r1%*%t(r1))
	  tempresult=optim(start_value,opt)
	  sol=tempresult$par
	  error<-tempresult$value
  } else
	{
	  start_value=rep(0,length(index.vars))
	  tempresult=optim(start_value,opt)
	  sol=tempresult$par
	  error<-tempresult$value
	}
  return(list(sol=sol,error=error))
}

