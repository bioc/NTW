`method.sse` <-
function(index.vars,X,pert)
{
  r1=X[index.vars,]
  rcd=1/kappa(r1%*%t(r1))
  if(rcd>=10^(-3))
  {
    sol = -pert %*% t(r1) %*% solve(r1%*%t(r1))
		error<-sum((-pert-sol%*%r1)^2)
  } else
  {
    sol = rep(0,length(index.vars))
		error <- Inf
  }
  return(list(sol=sol,error=error))
}

