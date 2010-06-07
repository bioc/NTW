\name{patterns.priorA}
\alias{backward}
\alias{forward}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Approaches to use a priori known gene association information for a single row estimation in matrix A }
\description{
  Given some information of the gene interaction network, for example the estimated gene association matrix \emph{pred.net} by other algorithms or from literature, NTW can use this prior information to enhance acurracy. NTW offers two approaches to infer the gene network, i.e. \emph{forward} and \emph{backward}, on the base of \emph{pred.net}. The former computes further  edges than the ones in \emph{pred.net}, while the latter prunes edges.
}
\usage{
backward(r, X, pert, topD, restk, cFlag, TA, noiseLevel)
forward(r, X, pert, topD, restk, cFlag, TA, noiseLevel)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{ A number to indicate the row of \emph{A} to be estimated when row \emph{r} of P is fixed. }
  \item{X}{ Gene expression data, a matrix with genes as rows and perturbations as columns. }
  \item{pert}{ Row \emph{r} in \emph{P}. }
  \item{topD}{ A parameter in NTW algorithm for keeping the top \emph{topD} combinations of non-zero regressors of row \emph{r} in \emph{A}, see \emph{vignette} for details. }
  \item{restk}{ The number of non-zero regressors for the estimation of row \emph{r} in \emph{A}. }
  \item{cFlag}{ A flag to identify the regression methods, "geo" for geometric mean method, "sse" for sum of square method and "ml" for maximum likelihood method. }
  \item{TA}{ A vector including the indexes of non-zero elements in row \emph{r} of the network containing a priori information, \emph{pred.net}.}
  \item{noiseLevel}{ Only used in "ml" method, to indicate the noise level in each perturbed experiment. }
}

\value{
   \item{ A.row }{ A vector of estimated row \emph{r}  in \emph{A}. }
  \item{ CrtValue }{ The minimum value from the objective function. }
  }
\author{ Wei Xiao, Yin Jin, Darong Lai, Xinyi Yang,Yuanhua Liu, Christine Nardini }

\keyword{ multivariate }