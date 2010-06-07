\name{AP.estimation.Srow}
\alias{AP.estimation.Srow}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Estimation of a single row in gene interaction matrix A and perturbation targets matrix P }
\description{
 Estimating a single row of gene interaction matrix \emph{A} and identifying the perturbations which target the corresponding gene of the row. For perturbation identifications, multiple perturbations are considered for one target (gene). Combinations of perturbations are first  assumed fixed. The single row in \emph{A} are then regressed according to equation AX=P with one of the three regression methods, \emph{geo} ,\emph{sse} and \emph{ml}. All these combinations will finally be optimized according to the difference between the predicted \emph{X} with the estimated \emph{A} and \emph{P} and the experimental values.
}
\usage{
AP.estimation.Srow(r, cMM.corrected, pred.net, X, IX, topD, restK, cFlag, sup.drop, numP, noiseLevel)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{ A number indicating the row of \emph{A} and \emph{P} to be estimated.  }
  \item{cMM.corrected}{ A flag to indicate whether a prior network is applied. }
  \item{pred.net}{ A matrix with  the same dimensions of \emph{A} for the prior network, which should be specified if cMM.corrected is 1, default is NULL. }
  \item{X}{ Gene expression data, a matrix with genes as rows and perturbations as columns. }
  \item{IX}{ The pre-estimated \emph{P} matrix according to the gene expression data \emph{X}, with the same dimensions of \emph{X} or \emph{P}. }
  \item{topD}{ A parameter in NTW algorithm for keeping the top topD combinations of non-zero regressors of row \emph{r} in \emph{A}, see \emph{vignette} for details. }
  \item{restK}{  A vector (length equals to \emph{nrow(A)}) with each element to indicate the number of non-zero regressors in the corresponding row of \emph{A}. }
  \item{cFlag}{ A flag to identify the regression methods, "geo" for geometric mean method, "sse" for sum of square method and "ml" for maximum likelihood method. }
  \item{sup.drop}{ A flag to show the pattern for using the prior gene association information.  \emph{1} for "forward" pattern and \emph{-1} for "backward" pattern. }
  \item{numP}{ A number set to limit the possibilities that one gene will be targeted by perturbations. That is at most \emph{numP} perturbations can directly perturb one gene. }
  \item{noiseLevel}{ Only used in "ml" method, to indicate the noise level in each perturbed experiment. }
}
 
\value{
  \item{ A.row }{ Estimation of the row \emph{r} in \emph{A}.}
  \item{ P.index }{ A vector to show which perturbations target the gene corresponding to row \emph{r}. }
 }

\author{ Wei Xiao, Yin Jin, Darong Lai, Xinyi Yang,Yuanhua Liu, Christine Nardini  }

\examples{
##single row estimation without prior gene association information, regression is done by "sse"##
data(sos.data)
X<-sos.data
X<-as.matrix(X)
IX<-P.preestimation(X, topK= round(2*nrow(X)))
restK=rep(ncol(X)-1, nrow(X))
topD = round(0.6*nrow(X))
topK = round(0.5*nrow(X))
numP = round(0.25*nrow(X))
result<-AP.estimation.Srow(r=1,cMM.corrected = 0, pred.net= NULL,X, IX,topD, restK,
               cFlag="sse",sup.drop = -1, numP, noiseLevel=0.1)
result$A.row
result$P.index

###single row estimation with prior gene association information, regression is done by "geo"###
pred.net<-matrix(round(runif(nrow(X)*nrow(X), min=0, max=1)), nrow(X), ncol(X))
result<-AP.estimation.Srow(r=1,cMM.corrected = 1, pred.net,X, IX,topD, restK,
             cFlag="geo",sup.drop = -1, numP, noiseLevel=0.1)
result$A.row
result$P.index 
}
\keyword{ arith }
