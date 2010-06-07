\name{NTW-package}
\alias{NTW-package}
\docType{package}
\title{ Gene interaction network and perturbation targets predictions }
\description{
This package includes the functions for estimating the gene-gene interaction network (a matrix, named \emph{A}, with genes as rows and columns) and the associated transcriptional targets of the perturbations (a matrix, named \emph{P}, with genes as rows and perturbations as columns). These estimations are computed with the NTW algorithm, a gene network inference algorithm based on ODE (ordinary differential equation) method, see \emph{reference}.  In this package, the whole \emph{A} matrix and \emph{P} matrix are estimated row by row with the function \emph{AP.estimation.Srow}, and built together with the function \emph{NTW}.  \emph{AP.estimation.Srow} can be used independently so that estimation of each row can be performed in parallel, improving computation time. For solving the steady state ODE equations, 3 regression methods are supplied: \emph{geo}, \emph{sse} and \emph{ml}, see details in the the corresponding function help pages.  In addition, in order to accelerate the estimation of matrix \emph{A}, an option is available to make use of some prior information such as gene association (output from other gene netwrok inference algorithms, or from literature) in \emph{NTW}. The regression methods used in forward or backward mode makes 6 possibilities available for estimating a single row of \emph{A} matrix.  The main functions in this package are listed below,
\itemize{
\item \emph{NTW}, to estimate the whole matrix \emph{A} and \emph{P} (if P is unknown).
\item \emph{AP.estimation.Srow}, to estimate  one single row in \emph{A}  and \emph{P} .
\item \emph{A.estimation.Srow}, to estimate  one single row in \emph{A}  with \emph{P} known.
\item \emph{backward} and \emph{forward}, to estimate one single row of matrix \emph{A} with different patterns of using prior gene association information. 
\item \emph{method.geo}, \emph{method.sse} and \emph{method.ml}, to estimate one single row of matrix \emph{A} with different regression methods.
\item \emph{comb.matrix}, sub-function to create all the combinations for regressor locations.
\item \emph{P.preestimation}, pre-estimate \emph{P} matrix according to the gene expression data.
}
}


\details{
\tabular{ll}{
Package: \tab NTW\cr
Type: \tab Package\cr
Version: \tab 0.99.0\cr
Date: \tab 2010-5-11\cr
License: \tab GPL-2 \cr 
LazyLoad: \tab yes\cr
}

}
\author{

Wei Xiao, Yin Jin, Darong Lai, Xinyi Yang, Yuanhua Liu, Christine Nardini

Maintainer: Yuanhua Liu <liuyuanhua@picb.ac.cn>
}

\references{
Applied method for the inference of gene networks: the bifidobacterium case.  to be submitted
}
\keyword{package}