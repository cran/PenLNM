\name{PenLNM.lmax}
\alias{PenLNM.lmax}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find the maximum penalty in the group l1 penalized LNM model
}
\description{
The function returns the maximum penalty value, which allows for only the intercept in the penalized LNM model. In estimating the LNM model, the implemented tuning parameter in the algorithm depends on the scale of variance-covariance matrix of underlying transformed compositions. Therefore, the maximum penalty could not be exactly identified before the model being estimated. The returned value is on the basis of the initial estimation of the variance-covariance matrix but is effective in practice. Starting at the returned maximum penalty value, the penalty grid points could be determined efficiently.
}

\usage{
PenLNM.lmax(base, X, W)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{base}{
The index of composition which will be used as the common denominator (base) in the log-ratio transformations 
}
  \item{X}{
The data matrix of covariates.
}
  \item{W}{
The data matrix of observed counts.
}
}

\value{
The empirical maximum penalty allows for only the intercept in the LNM model
}

\references{
Fan Xia, Jun Chen, Wing Kam Fung, Hongzhe, Li (2012). A Logistic Normal Multinomial Regression Model for Microbiome Compositional Data Analysis. (Submitted)
}
\author{
FAN XIA <phoebexia@yahoo.com>
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{PenLNM.CV}}
}
\examples{
data(Nutrients)
data(Genera3)

lamda.max <- PenLNM.lmax(2,Nutrients,Genera3)

## Start at lamda.max, a grid of penalties are determined, e.g.,

lamda.set <- lamda.max*0.9^(0:9)

## for 10 grid points

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ tuning }
\keyword{ penalty }% __ONLY ONE__ keyword per line