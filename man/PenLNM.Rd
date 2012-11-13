\name{PenLNM}
\alias{PenLNM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate the group l1 penalized LN/LNM model.
}
\description{
This is the major function in the package. It estimates the penalized models and return parameter estimates.
}
\usage{
PenLNM(base, X, W, lamda, model = c("LN", "LNM"), Iter1 = 10, 
burn = 5, Iter2 = 1000, MH.burn = 500, 
MH.scale = 0.5, save = TRUE, plot = TRUE)
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
The data matrix of observed counts.
}
  \item{lamda}{
The selected penalty
}
  \item{model}{
The method used to model the compositional data set. Default is the simple "LN" model, alternatively, model="LNM" is specified.
}
  \item{Iter1}{
The total number of EM iterations.
}
  \item{burn}{
The number of EM iterations for burn-in.
}
  \item{Iter2}{
The total number of Metropolis-Hasting iterations in each EM iteration, default value is 1000.
}
  \item{MH.burn}{
The number of Metropolis-hasting samples for burn in, default value is 500
}
  \item{MH.scale}{
The step size in the Metropolis Hasting algorithm and the default value is 0.5. The value of MH.scale will affect the acceptance ratio in MH algorithm. Users are suggested to set a MH.scale that generates averaged acceptance ratio around 0.3-0.5. Empirically, MH.scale ranging from 0.5-1 is found to be effective in generating efficient MH steps.
}
  \item{save}{
Whether or not to save the intermediate results in the implementation of the MCEM algorithm. If save = TRUE (default), A LIST containing the estimates in each EM iteration of regression constant vector, regression coefficient matrix, variance-covariance matrix of the underlying transformed compositions, and the averaged acceptance ratios of the MH steps  will be stored in a .RData file named "listfull.RData"; otherwise, no intermediate results will be saved.
}
  \item{plot}{
Whether or not tp show the fitting plot, default setting is plot = TRUE.
}
}

\value{
Return a LIST containing:
\item{lamda}{The penalty parameter value.}
\item{master_pred}{The index of selected covariates.}
\item{b0.simplex}{The vector of logistic transformed LN/LNM regression constants, which could be seen as the baseline composition.}
\item{b.simplex}{The matrix of logistic transformed LN/LNM regression coefficients, which could be seen as the perturbations from selected covariates.}
\item{MSPE}{The mean square prediction error based on fitted counts and the observed counts.}
}
\references{
Fan Xia, Jun Chen, Wing Kam Fung, Hongzhe, Li (2012). A Logistic Normal Multinomial Regression Model for Microbiome Compositional Data Analysis. (Submitted)
}
\author{
FAN XIA <phoebexia@yahoo.com>
}

\seealso{
\code{\link{PenLNM.CV}}
}

\examples{
data(Nutrients)
data(Genera3)

\dontrun{
PenLNM(base=2,Nutrients,Genera3,lamda=35,model="LNM")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ logistic-normal-multinomial }
\keyword{ regression }
\keyword{ variable-selection }% __ONLY ONE__ keyword per line
