\name{PenLNM.CV}
\alias{PenLNM.CV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the cross-validation of the penalized LNM approach for a pre-assigned penalty value.
}
\description{
In practice, the calculation of cross-validation for the LNM model is time-consuming that parallel programing is strongly suggested. 
}
\usage{
PenLNM.CV(base, X, W, lamda, Iter1 = 10, burn = 5, Iter2 = 1000, 
MH.burn = 500, MH.scale = 0.5, cv.folds = 5)
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
The pre-assigned penalty
}
  \item{Iter1}{
The total number of EM iterations.
}
  \item{burn}{
The number of EM iterations for burn-in.
}
  \item{Iter2}{
The total number of Metropolis-hasting iterations in each EM iteration, default is 1000.
}
  \item{MH.burn}{
The number of Metropolis-hasting samples for burn in, default value is 500
}
  \item{MH.scale}{
The step size in the Metropolis Hasting algorithm, default value is 0.5. The value of MH.scale will affect the acceptance ratio in MH algorithm. Users are suggested to set a MH.scale that generates averaged acceptance ratio around 0.3-0.5. Empirically, MH.scale ranging from 0.5-1 is found to be effective in generating efficient MH steps.
}
  \item{cv.folds}{
The number of folds in the calculation of cross validation, default value is 5.
}
}

\value{
Return a vector with the first element being the calculated cross-validation and the second element being the corresponding penalty. 
%% ...
}
\references{
Fan Xia, Jun Chen, Wing Kam Fung, Hongzhe, Li (2012). A Logistic Normal Multinomial Regression Model for Microbiome Compositional Data Analysis. (Submitted)
}
\author{
FAN XIA <phoebexia@yahoo.com>
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{PenLNM.lmax}}
}
\examples{

data(Nutrients)
data(Genera3)

N.lamda <- 10 
lamda.max <- PenLNM.lmax(2,Nutrients,Genera3)
lamda.set <- lamda.max*0.9^(0:(N.lamda-1))

\dontrun{
# The package "Rmpi" facility the parallel programing in R.
# Each slave estimate the cross-validation of one penalty in the penalty grid.

## Load the R MPI package if it is not already loaded.
if (!is.loaded("mpi_initialize")) {
 library("Rmpi")
 }

## Spawn N.lamda slaves 
mpi.spawn.Rslaves(N.lamda)

# In case R exits unexpectedly, have it automatically clean up
# resources taken up by Rmpi (slaves, memory, etc...)
.Last <- function(){
  if (is.loaded("mpi_initialize")){
    if (mpi.comm.size(1) > 0){
      print("Please use mpi.close.Rslaves() to close slaves.")
      mpi.close.Rslaves()
    }
    print("Please use mpi.quit() to quit R")
    .Call("mpi_finalize")
  }
}

foldslave=function(){
  library(remMap)
  library(MASS)
  library(Matrix)
  library(PenLNM)
  
  data(Nutrients)
  data(Genera3)

  lamda <- lamda.set[foldNumber]
  PenLNM.CV(base=2,Nutrients,Genera3,lamda)
}

## send the commands to slaves
mpi.bcast.Robj2slave(lamda.set)
mpi.bcast.cmd(foldNumber <- mpi.comm.rank())
mpi.bcast.Robj2slave(foldslave)

Res <- mpi.remote.exec(foldslave())
print(Res)

## Based on Res, select the lamda with the corresponding smallest cross validation.

mpi.close.Rslaves()
mpi.quit(save = "no")
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ parallel }
\keyword{ cross-validation }% __ONLY ONE__ keyword per line
