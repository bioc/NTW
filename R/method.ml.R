`method.ml` <-
function(index.vars,X,pert,noiseLevel)
{
	opt=function(bta)
	{
	  errorY<-(-pert-t(bta)%*%r1)
		variancematrix<-diag(as.vector((t(r1))^2%*%(bta)^2)+1)
		return(1/dmvnorm(errorY,rep(0,length(errorY)),((noiseLevel)^2)*variancematrix))
  }
	r1=X[index.vars,]
	rcd=1/kappa(r1%*%t(r1))
	if(rcd>=10^(-9))
	{
	  start_value=t(-pert %*% t(r1) %*% solve(r1%*%t(r1)))
	  tempresult=optim(start_value,opt)
	  sol=tempresult$par
	  error<-tempresult$value
	} else
	{
    sol<-rep(0,length(index.vars))
	  error<-Inf
	}
  return(list(sol=sol,error=error))
}

