\name{P.preestimation}
\alias{P.preestimation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Pre-estimation of the transcriptional pertubation targets matrix P }
\description{
  Pre-estimate the potential transcriptional perturbation targets matrix \emph{P} according to gene expression data \emph{X} . Those genes with the changes in top topK will be assumed as possible targets of the perturbations. 
}
\usage{
P.preestimation(X, topK)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{ Gene expression data, a matrix with rows as genes and columns as experiments. }
  \item{topK}{ The number of possible targets of the perturbations. }
}
\value{
A matrix with the same structure of \emph{X} or \emph{P}.
}
\author{ Wei Xiao, Yin Jin, Darong Lai, Xinyi Yang, Yuanhua Liu, Christine Nardini }
\examples{
data(sos.data)
X<-sos.data
X<-as.matrix(X)
IX<-P.preestimation(X, topK= round(0.6*nrow(X)))
IX
}
\keyword{ arith }