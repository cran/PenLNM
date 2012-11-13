\name{PenLNM-package}
\alias{PenLNM-package}
\docType{package}
\title{
Group $l_{1}$ penalized logistic normal multinomial \strong{LNM} regression model.
}
\description{
The package is developed for the cross-sectional study of diet and stool microbiome data. The model facility variable selection and therefore accommodate the nutrients of high-dimensionality. The group $l_{1}$ penalized \strong{LN} approach is also implemented in the package. In addition, the package contains functions for tuning the penalty parameter. 
}

\author{
FAN XIA <phoebexia@yahoo.com>
}
\references{
Fan Xia, Jun Chen, Wing Kam Fung, Hongzhe, Li (2012). A Logistic Normal Multinomial Regression Model for Microbiome Compositional Data Analysis. (Submitted)
}

\keyword{ package }
\keyword{ compositional-data }
\keyword{ logistic-normal-multinomial }
\keyword{ regression }
\keyword{ variable-selection }

\examples{
data(Genera3)
data(Nutrients)

# Analyze the data with group l1 penalized LN model
# Tuning parameter in the LN model
lamda.LN <- PenLN.CV(base=2,Nutrients,Genera3)[[2]]

# Fit the model with selected penalty
PenLNM(base=2,Nutrients,Genera3,lamda=lamda.LN,model="LN")

# Analyze the data with group l1 penalized LNM model

# Find the maximum tuning parameter
lamda.max <- PenLNM.lmax(2,Nutrients,Genera3)

# Determine a grid of penalties for tuning starting at lamda.max, e.g.,
lamda.set <- lamda.max*0.9^(0:9)

# Since the estimation of LNM model is time consuming, parallel programing is strongly suggested for tuning penalty parameter in the LNM model
# See details in the page of function PenLNM.CV()

# Fit the model with selected penalty, e.g., 35.
\dontrun{
PenLNM(base=2,Nutrients,Genera3,lamda=35,model="LNM")
}

#It takes around 30 minutes to 1 hour to generate results for the default settings with 20 EM iterations and 1000 MH iterations in each EM.
}


