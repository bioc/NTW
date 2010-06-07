\name{NTW}
\alias{NTW}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Estimation of gene interaction matrix A and perturbation targets matrix P }
\description{
 This function is used to estimate the whole gene interaction matrix \emph{A} and the perturbation targets matrix \emph{P}, row-wise, using the NTW algorithm (see \emph{reference}), based on ODE method.  In this method, the linearized ODE can solved using 3 regression methods: \emph{geo}, \emph{sse} and \emph{ml}. In order to save computation time, and improve results, NTW offers the opportunity to input gene association information output from other algorithms or from the literature. The non-null regressors in the gene association network will help fix the regressors to be estimated in the final matrix \emph{A}. Two ways are supplied to use the non-zero information, namely \emph{forward} and \emph{backward} approaches. In the "backward" pattern, only the non-zero positions in the prior gene association network will be used as regressors in \emph{A}. While in the "forward" pattern, both these non-zero positions and some other possible positions (depending on \emph{restK} ) in \emph{A} are used as regressors.
}
\usage{
NTW(X, restK, topD, topK = NULL, P.known = NULL, cFlag, pred.net = NULL, sup.drop = -1, numP = NULL, noiseLevel = 0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{ Gene expression data, a matrix with genes as rows and perturbations as columns. }
  \item{restK}{ A vector (length equals to nrow(A)) with each element to indicate the number of non-zero regressors in the corresponding row of \emph{A}. }
  \item{topD}{ A parameter in NTW algorithm for keeping the top \emph{topD} combinations of non-zero regressors of a single row in \emph{A}, see \emph{vignette} for details. }
  \item{topK}{ The number of possible targets of the perturbations, used for pre-estimate the perturbation targets matrix \emph{P} . }
  \item{P.known}{ A known P matrix with the same dimensions of \emph{X}.}
  \item{cFlag}{ A flag to tell the regression methods, "geo" for geometric mean method, "sse" for sum of square method and "ml" for maximum likelihood method. }
  \item{pred.net}{ A matrix with the same dimensions of \emph{A} for the prior gene association information. Default is NULL. }
  \item{sup.drop}{ An indication to show the pattern for using the prior gene association information.  \emph{1} for "forward" pattern and \emph{-1} for "backward" pattern. }
  \item{numP}{ A number set to limit the possibilities that one gene will be targeted by perturbations. That is at most \emph{numP} perturbations can directly perturb one gene. }
  \item{noiseLevel}{ Only used in \emph{ml} method, to indicate the noise level in each perturbed experiment. }
}
 
\value{
  \item{ est.A }{ Estimated gene interaction matrix \emph{A}, with genes as rows and columns.}
  \item{ est.P }{ Estimated perturbation targets matrix \emph{P}, with genes as rows and perturbations as columns. }
}
\references{
Applied method for the inference of gene networks: the bifidobacterium case.  to be submitted
}
\author{ Wei Xiao, Yin Jin, Darong Lai, Xinyi Yang,Yuanhua Liu, Christine Nardini }
\examples{

##NTW testing without prior gene association information, regression is done by "sse"##
data(sos.data)
X<-sos.data
X<-as.matrix(X)
restK=rep(ncol(X)-1, nrow(X))
topD = round(0.6*nrow(X))
topK = round(0.5*nrow(X))
numP = round(0.25*nrow(X))
result<-NTW(X, restK, topD, topK, P.known=NULL, cFlag="sse", 
              pred.net = NULL, sup.drop = -1,numP, noiseLevel=0.1)
result$est.A
result$est.P

##NTW testing with prior gene association information, regression is done by "geo"##
pred.net<-matrix(round(runif(nrow(X)*nrow(X), min=0, max=1)), nrow(X), nrow(X))
result<-NTW(X, restK, topD, topK, P.known=NULL, cFlag="geo",
              pred.net, sup.drop = -1,numP, noiseLevel=0.1)
result$est.A
result$est.P

}
\keyword{ arith }
