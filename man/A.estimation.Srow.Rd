\name{A.estimation.Srow}
\alias{A.estimation.Srow}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Estimation of a single row in matrix A with the perturbation targets matrix P known }
\description{
Estimating a single row of gene interaction matrix \emph{A} when the perturbation targets matrix P is given. The single row in \emph{A} is then regressed according to the equation AX=P with one of the three regression methods, \emph{geo}, \emph{sse} and \emph{ml} .
}
\usage{
A.estimation.Srow(r, cMM.corrected, pred.net, X, P.known, topD, restK, cFlag, sup.drop, noiseLevel)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{ A number indicating the row of \emph{A} to be estimated. }
  \item{cMM.corrected}{ A flag to indicate whether a prior network is applied. }
  \item{pred.net}{ A matrix with  the same dimensions of \emph{A} for prior network, which should be specified if cMM.corrected is 1, default is NULL. }
  \item{X}{ Gene expression data, a matrix with genes as rows and perturbations as columns. }
  \item{P.known}{ A known P matrix with the same dimensions of \emph{X}. }
  \item{topD}{ A parameter in NTW algorithm for keeping the top \emph{topD} combinations of non-zero regressors of row \emph{r} in \emph{A}, see \emph{vignette} for details. }
  \item{restK}{ A vector (length equals to \emph{nrow(A)}) with each element to indicate the number of non-zero regressors in the corresponding row of \emph{A}. }
  \item{cFlag}{ A flag to tell the regression methods, "geo" for geometric mean method, "sse" for sum of square method and "ml" for maximum likelihood method. }
  \item{sup.drop}{ An indication to identify the pattern for using the prior gene association information.  \emph{1} for "forward" pattern and \emph{-1} for "backward" pattern, see \emph{vignette} for details. }
  \item{noiseLevel}{ Only used in "ml" method, to indicate the noise level in each perturbed experiment.}
}

\value{
 \item{ A.row }{ A vector of estimated row \emph{r} in \emph{A}.}
}
\author{ Wei Xiao, Yin Jin, Darong Lai, Xinyi Yang,Yuanhua Liu, Christine Nardini }

\examples{
##single row estimation without prior gene association information, regression is done by "sse"##
data(sos.data)
X<-sos.data
X<-as.matrix(X)
P.known<-matrix(round(runif(nrow(X)*ncol(X), min=0, max=1)), nrow(X), ncol(X))
restK=rep(ncol(X)-1, nrow(X))
topD = round(0.6*nrow(X))
topK = round(0.5*nrow(X))
result<-A.estimation.Srow(r=1,cMM.corrected = 0, pred.net= NULL,X,P.known, topD, restK, 
              cFlag="sse",sup.drop = -1, noiseLevel=0.1)
result$A.row

##single row estimation with prior gene association information, regression is done by "geo"###
pred.net<-matrix(round(runif(nrow(X)*nrow(X), min=0, max=1)), nrow(X), ncol(X))
result<-A.estimation.Srow(r=1,cMM.corrected = 1, pred.net,X,P.known,topD, restK,
             cFlag="geo",sup.drop = -1, noiseLevel=0.1)
result$A.row
}
\keyword{ arith }