\name{PenLN.CV}
\alias{PenLN.CV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Perform penalized LN model over a grid of penalty parameter values and returns the selected penalty based on cross-validation criterion.
}
\description{
The cross validation is used to determined the penalty parameter value in the naive LN approach. Compared to the estimation of penalized LNM model where the MCEM algorithm is implemented, the estimation of penalized LN model costs little computation efforts.
}

\usage{
PenLN.CV(base, X, W, N.lamda = 20, lamda.distance = 0.9, cv.folds = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{base}{
The index of taxa which will be used as the common denominator (base) in the log-ratio transformations
}
  \item{X}{
The data matrix of covariates.
}
  \item{W}{
The data matrix of taxa counts.
}
  \item{N.lamda}{
The number of penalty grid points
}
  \item{lamda.distance}{
Start at the maximum penalty, the distances between different grid points. The maximum penalty value will be computed automatically and then the grid of penalties are given by: lamda.max*lamda.distance^(0:(N.lamda-1)).
}
  \item{cv.folds}{
The number of folds in the calculation of cross validation, default value is 5.
}
}

\value{
Return a LIST containing:
\item{"Cross validation of different lamda"}{The vector of cross validations for different penalties}
\item{"the selected lamda"}{The selected penalty}
}

\references{
Fan Xia, Jun Chen, Wing Kam Fung, Hongzhe, Li (2012). A Logistic Normal Multinomial Regression Model for Microbiome Compositional Data Analysis. (Submitted).
}
\author{
FAN XIA <phoebexia@yahoo.com>
}

\examples{
data(Nutrients)
data(Genera3)

lamda.LN <- PenLN.CV(base=2,Nutrients,Genera3)[[2]]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ logistic-normal }
\keyword{ cross-validation }% __ONLY ONE__ keyword per line
