\name{methods.regression}
\alias{method.geo}
\alias{method.sse}
\alias{method.ml}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Regression methods to estimate a single row in A with fixed perturbations}
\description{
  Regression methods to estimate a single row in a gene interaction network (\emph{A}) with perturbations (\emph{P}) fixed. These methods differ at the regression criterions: the objective function of \emph{geo} is the geometric mean, sum of square for \emph{sse} and maximum likelihood for \emph{ml}, see \emph{vignette} for the details.
}

\usage{
method.geo(index.vars, X, pert)
method.sse(index.vars, X, pert)
method.ml(index.vars, X, pert, noiseLevel)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{index.vars}{ A vector to select the rows in \emph{X} on which regression will be used. It sets a group of combination of rows in \emph{X}. }
  \item{X}{ Gene expression data, a matrix with genes as rows and perturbations as columns. }
  \item{pert}{ A vector of row \emph{r} in \emph{P}.}
   \item{noiseLevel}{ Indicate the noise level in each perturbed experiment.}
  }

\value{
   \item{ sol }{ A vector of regression result }
  \item{ error }{ The result of the objective function. }
}
\author{Wei Xiao, Yin Jin, Darong Lai, Xinyi Yang,Yuanhua Liu, Christine Nardini}

\keyword{ multivariate }
