\name{comb.matrix}
\alias{comb.matrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Create all combinations of vectors }
\description{
  Create all combinations of vectors especially for matrices. 
}
\usage{
comb.matrix(x, y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A vector. }
  \item{y}{ A vector. }
}
\value{
  A matrix  with  nrow(x)*nrow(y) rows and ncol(x)+ncol(y) columns.
}
\author{ Wei Xiao, Yin Jin, Darong Lai, Xinyi Yang, Yuanhua Liu, Christine Nardini }
\examples{
###A matrix with only one row is obtained####
x<-c(1,2,3)
y<-c(4,5)
comb.matrix(x,y)
###A matrix with 2 rows and 4 columns is obtained####
x<-matrix(x,1,)
y<-matrix(y,,1)
comb.matrix(x,y)

}
\keyword{ arith }

