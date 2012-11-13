########################## PenLNM ################################
## PenLNM: Group l1 penalized LNM model for cross-sectional study 
## of diet and stool microbiome composition.

## Author: FAN XIA <phoebexia@yahoo.com

## Date: 14 Oct 2012 

## Reference: Fan Xia, Jun Chen, Wing Kam Fung, Hongzhe, Li (2012). 
## A Logistic Normal Multinomial Regression Model for Microbiome
## Compositional Data Analysis. (Submitted)

########################## PenLNM ################################
## require the following packages
require(remMap)
require(MASS)
require(Matrix)
###################################################################
PenLN.CV <- function(base,X,W,N.lamda=20,lamda.distance=0.9,cv.folds=5){
  ## function to tuning the penalty in the penalized LN model
  ## N.lamda: number of grid points
  ## lamda.distance: distance between grid points
  Ini.Logratios <- function(base,W){
    ## functions for initialize the logratio responses in LN distribution
    ## base: the base taxa used in logratio transformation
    W.m <- as.matrix(W)
    N <- nrow(W.m)
    Q <- ncol(W.m) #dim of covariates and count vectors
    M=apply(W.m,1,sum)
    
    Psedo_Z.m <- matrix(0,N,Q) 
    Y.m <- matrix(0,N,Q-1) 
    
    for(i in 1:N){
      zeros <- which(W.m[i,]==0)
      nzeros <- which(W.m[i,]!=0)
      Psedo_Z.m[i,zeros] <- (W.m[i,zeros]+0.05)/(M[i]+0.05*length(zeros))
      Psedo_Z.m[i,nzeros]<- (W.m[i,nzeros])/(M[i]+0.05*length(zeros))
      Y.m[i,] <- log(Psedo_Z.m[i,-base]/Psedo_Z.m[i,base])
    }
    attr(Y.m,"center")=apply(Y.m,2,mean)
    attr(Y.m,"base")=paste("the",base,"th taxa")
    Y.m
  }
  
  X.m <- scale(X);W.m <- as.matrix(W)
  Y.m <- Ini.Logratios(base,W)
  
  N <- nrow(W.m)
  P <- ncol(X.m)
  Q <- ncol(W.m)
  M <- apply(W.m,1,sum)
  
  lam_maxi <- NULL
  CY.m <- scale(Y.m,scale=FALSE)
  for(i in 1:P){
    lam_maxi[i] <- sqrt(sum((t(X.m[,i])%*%CY.m)^2))
  }
  lam_max <- max(lam_maxi)
  lamda.vector <- lam_max*lamda.distance^(0:(N.lamda-1))
  
  MSPE_CV=matrix(0,N.lamda,cv.folds)
  ##now, specify the traning and testing data set and seperately center the 
  ##covariates and responses
  for(la in 1:N.lamda){
    lamda <- lamda.vector[la]
    
    for(cv in 1:cv.folds){
      test.ID <- floor(N*(cv-1)/cv.folds+1):floor(N*cv/cv.folds)
      M_test <- M[test.ID]
      W.m_test <- W.m[test.ID,]
      Y.m_test <- scale(Y.m[test.ID,],scale=FALSE)
      X.m_test <- scale(X.m[test.ID,])
      N_test <- length(test.ID)
      
      M_train <- M[-test.ID]
      W.m_train <- W.m[-test.ID,]
      Y.m_train <- scale(Y.m[-test.ID,],scale=FALSE)
      X.m_train <- scale(X.m[-test.ID,])
      N_train <- N-N_test
      
      #fit the model with training dataset for a specific lamda
      b_ini <- remMap(X.m_train,Y.m_train,0,lamda)$phi
      master_pred <- which(apply(abs(b_ini),1,sum)!=0)
      
      if(length(master_pred)==0){
        eY.m <- rep(1,N_test)%*%t(attr(Y.m_test,"scaled:center"))
      }
      
      if(length(master_pred)>0){
        MLE_fit <- lm(Y.m_train~X.m_train[,master_pred])
        b <- MLE_fit$coefficients[-1,]
        eY.m <- as.matrix(X.m_test[,master_pred])%*%b+rep(1,N_test)%*%t(attr(Y.m_test,"scaled:center"))
      }
      exp_eY <- exp(eY.m)
      sum_exp_eY <- apply(exp_eY,1,sum) 
      e_Z.m <- matrix(0,N_test,Q) 
      e_W.m <- matrix(0,N_test,Q)
      for(i in 1:N_test){
        e_Z.m[i,-base] <- exp_eY[i,]/(sum_exp_eY[i]+1)
        e_Z.m[i,base] <- 1/(sum_exp_eY[i]+1)
        e_W.m[i,] <- M_test[i]*e_Z.m[i,]
      }
      MSPE_CV[la,cv] <- sqrt(sum((e_W.m-W.m_test)^2))
    }
  }
  
  cv.all <- apply(MSPE_CV,1,mean)
  names(cv.all) <- lamda.vector
  lamda.final <- lamda.vector[which.min(cv.all)]
  Result <- list(cv.all,lamda.final)
  names(Result) <- c("Cross validation of different lamda",
                       "the selected lamda")
  Result
}

PenLNM.lmax <- function(base,X,W){
  ## Find the maximum value of penalty for the LNM model
  Ini.Logratios <- function(base,W){
    ## functions for initialize the logratio responses in LN distribution
    ## base: the base taxa used in logratio transformation
    W.m <- as.matrix(W)
    N <- nrow(W.m)
    Q <- ncol(W.m) #dim of covariates and count vectors
    M=apply(W.m,1,sum)
    
    Psedo_Z.m <- matrix(0,N,Q) 
    Y.m <- matrix(0,N,Q-1) 
    
    for(i in 1:N){
      zeros <- which(W.m[i,]==0)
      nzeros <- which(W.m[i,]!=0)
      Psedo_Z.m[i,zeros] <- (W.m[i,zeros]+0.05)/(M[i]+0.05*length(zeros))
      Psedo_Z.m[i,nzeros]<- (W.m[i,nzeros])/(M[i]+0.05*length(zeros))
      Y.m[i,] <- log(Psedo_Z.m[i,-base]/Psedo_Z.m[i,base])
    }
    attr(Y.m,"center")=apply(Y.m,2,mean)
    attr(Y.m,"base")=paste("the",base,"th taxa")
    Y.m
  }
  
  X.m <- scale(X);W.m <- as.matrix(W)
  P <- ncol(X.m)
  Y.m <- Ini.Logratios(base,W)
  CY.m <- scale(Y.m,scale=FALSE)
  sigma <- var(CY.m)
  
  lam_maxi <- NULL
  U <- chol(solve(sigma))
  CYu.m <- as.matrix(scale(CY.m%*%U,center=TRUE,scale=FALSE)) 
  
  lam_maxi <- NULL
  for(i in 1:P){
    lam_maxi[i] <- sqrt(sum((t(X.m[,i])%*%CYu.m)^2))
  }
  lam_max <- max(lam_maxi)
  lam_max
}

## We require the parallel programing for tuning the LNM model 
## In R, the package "Rmpi" encourages parallel programing
## the grid points of penalty start at PenLNM.lmax 

PenLNM.CV <- function(base,X,W,lamda,Iter1=10,burn=5, Iter2=1000,
                          MH.burn=500,MH.scale=0.5,cv.folds=5)
{
  Ini.Logratios <- function(base,W){
    ## functions for initialize the logratio responses in LN distribution
    ## base: the base taxa used in logratio transformation
    W.m <- as.matrix(W)
    N <- nrow(W.m)
    Q <- ncol(W.m) #dim of covariates and count vectors
    M=apply(W.m,1,sum)
    
    Psedo_Z.m <- matrix(0,N,Q) 
    Y.m <- matrix(0,N,Q-1) 
    
    for(i in 1:N){
      zeros <- which(W.m[i,]==0)
      nzeros <- which(W.m[i,]!=0)
      Psedo_Z.m[i,zeros] <- (W.m[i,zeros]+0.05)/(M[i]+0.05*length(zeros))
      Psedo_Z.m[i,nzeros]<- (W.m[i,nzeros])/(M[i]+0.05*length(zeros))
      Y.m[i,] <- log(Psedo_Z.m[i,-base]/Psedo_Z.m[i,base])
    }
    attr(Y.m,"center")=apply(Y.m,2,mean)
    attr(Y.m,"base")=paste("the",base,"th taxa")
    Y.m
  }
  
  OLS_fun <- function(X.m,Y.m,b_ini){
    ## function specially for re-estimation
    P <- ncol(X.m); Q <- ncol(Y.m)+1
    CY.m <- scale(Y.m,scale=FALSE)
    b <- matrix(0,P,Q-1)
    for(q in 1:(Q-1)){
      b_q <- b_ini[,q]
      edges <- which(b_q!=0) 
      #the estimated edges for yq
      if(length(edges)>0){
        b_q[edges] <- ginv(t(X.m[,edges])%*%X.m[,edges])%*%t(X.m[,edges])%*%CY.m[,q]
      }
      b[,q]<- b_q
    }
    b
  }
  
  X.m <- scale(X);W.m <- as.matrix(W)
  Y.m <- Ini.Logratios(base,W)
  
  N <- nrow(W.m)
  P <- ncol(X.m)
  Q <- ncol(W.m)
  M <- apply(W.m,1,sum)
  
  MSPE_CV <- NULL
  
  for(cv in 1:cv.folds){
    cat(paste("the fold is",cv,"\n"))
    test.ID <- floor(N*(cv-1)/cv.folds+1):floor(N*cv/cv.folds)
    M_test <- M[test.ID]
    W.m_test <- W.m[test.ID,]
    Y.m_test <- Y.m[test.ID,]
    CY.m_test <- scale(Y.m_test, scale=FALSE)
    X.m_test <- scale(X.m[test.ID,])
    N_test <- length(test.ID)
    
    M_train <- M[-test.ID]
    W.m_train <- W.m[-test.ID,]
    Y.m_train <- Y.m[-test.ID,]
    CY.m_train <- scale(Y.m_train, scale=FALSE)
    X.m_train <- scale(X.m[-test.ID,])
    N_train <- N-N_test
    
    ## initialize the parameter
    b0 <- attr(CY.m_train,"scaled:center")
    sigma <- var(Y.m_train)
    U <- chol(solve(sigma))
    CYu.m <- as.matrix(scale(CY.m_train%*%U,scale=FALSE)) 
    b_ini <- remMap(X.m_train,CYu.m,0,lamda)$phi 
    b <- OLS_fun(X.m_train,Y.m_train,b_ini)
    eY.m_train <- X.m_train %*% b+rep(1,N_train)%*%t(b0)
    
    list.b <- list(b)
    list.b0 <- list(b0)
    list.sigma <- list(sigma)
    Y.m_train_em <- Y.m_train
    
    ## Start the MCEM algorithm
    for(em in 1:Iter1){
      MH_samples <- list(NULL)
      ## store all the MH samples in one EM step
      for(i in 1:N_train){
        ## MH independently proceeded for each individual
        ## therefore (N*Iter2) MH iterations in total in each EM step
        YI.store <- matrix(0,Iter2,Q-1)
        
        YI <- Y.m_train_em[i,] 
        WI <- W.m_train[i,]
        muI <- eY.m_train[i,]
        
        for(mc in 1:Iter2){
          YI_n <- YI+rnorm(Q-1,0,MH.scale)
          
          Log_pt1 <- sum(WI)*(log(sum(exp(YI))+1)
                              -log(sum(exp(YI_n))+1))
          Log_pt2 <- sum(WI[-base]*(YI_n-YI))
          Log_pt3 <- -0.5*t(YI_n-muI)%*%solve(sigma)%*%(YI_n-muI)
          Log_pt4 <- -0.5*t(YI-muI)%*%solve(sigma)%*%(YI-muI)
          R_ratio <- Log_pt1+Log_pt2+Log_pt3-Log_pt4
          acceptance <- min(1,exp(R_ratio))
          ## the acceptance ratio
          temp<-runif(1)
          if(temp < acceptance){YI <- YI_n}
          YI.store[mc,]<- YI
        }
        Y.m_train_em[i,]<- apply(YI.store[(MH.burn+1):Iter2,],2,mean)
        ## update the Y.m based on the stablized MH samples
        MH_samples <- c(MH_samples,list(YI.store))
      }
      
      ## update b0 and sigma with the MH-samples
      SV_sum<-matrix(0,Q-1,Q-1)
      for(mc in (MH.burn+1):Iter2){
        Y.m_mc <- matrix(0,N_train,Q-1)
        for(i in 1:N_train){
          Y.m_mc[i,] <- MH_samples[[i+1]][mc,]
        }
        eps_mc <- Y.m_mc-eY.m_train
        SV_mc <- t(eps_mc)%*%eps_mc
        SV_sum <- SV_sum+SV_mc
      }
      
      sigma.new <- SV_sum/(N_train*(Iter2-MH.burn))
      b0.new <- apply(Y.m_train_em,2,mean)
      
      ##update b
      U <- chol(solve(sigma.new))
      CY.m <- Y.m_train_em-rep(1,N_train)%*%t(b0.new)
      CYu.m <- as.matrix(scale(CY.m%*%U,scale=FALSE)) 
      b_ini <- remMap(X.m_train,CYu.m,0,lamda)$phi
      
      b.new <- OLS_fun(X.m_train,Y.m_train_em,b_ini)
      
      sigma <- sigma.new
      b0 <- b0.new
      b <- b.new
      eY.m_train <- X.m_train%*%b+rep(1,N_train)%*%t(b0)
      
      list.b <- c(list.b,list(b))
      list.b0 <- c(list.b0,list(b0))
      list.sigma <- c(list.sigma,list(sigma))    
    }
    
    sigma.sum <- matrix(0,Q-1,Q-1)
    b0.sum <- rep(0,Q-1)
    
    for(em in (burn+1):Iter1){
      sigma.sum <- sigma.sum+list.sigma[[em+1]]
      b0.sum <- b0.sum+list.b0[[em+1]]
    }
    
    sigma.mean <- sigma.sum/(Iter1-burn)
    b0.mean <- b0.sum/(Iter1-burn)
    
    b.med <- matrix(0,P,Q-1)
    # b.mean=matrix(0,P,Q-1)
    for(i in 1:P){
      for(j in 1:(Q-1)){
        bij_store <- NULL
        for(em in (burn+1):Iter1){
          bij_store <- c(bij_store,list.b[[em+1]][i,j])
        }
        b.med[i,j] <- median(bij_store)
        # b.mean[i,j]=mean(b.mean)
      }
    }
    
    master_pred <- which(apply(abs(b.med),1,sum)!=0)
    
    if(length(master_pred)== 0){
      eY.m <- rep(1,N_test)%*%t(attr(CY.m_test,"scaled:center"))
    }
    
    if(length(master_pred)>0){
      eY.m <- as.matrix(X.m_test[,master_pred]) %*% b.med[master_pred,]+rep(1,N_test)%*%t(attr(CY.m_test,"scaled:center"))
    }
    
    exp_eY <- exp(eY.m)
    sum_exp_eY <- apply(exp_eY,1,sum) 
    e_Z.m <- matrix(0,N_test,Q) 
    e_W.m <- matrix(0,N_test,Q)
    for(i in 1:N_test){
      e_Z.m[i,-base] <- exp_eY[i,]/(sum_exp_eY[i]+1)
      e_Z.m[i,base] <- 1/(sum_exp_eY[i]+1)
      e_W.m[i,] <- M_test[i]*e_Z.m[i,]
    }
    MSPE_CV[cv] <- sqrt(sum((e_W.m-W.m_test)^2))
  }
  
  CV <- mean(MSPE_CV)
  Result <- c(CV,lamda)
  names(Result) <- c("Cross validation","lamda")
  Result
}


PenLNM <- function(base,X,W,lamda,model=c("LN","LNM"),
                   Iter1=10,burn=5,Iter2=1000,MH.burn=500,MH.scale=0.5,save=TRUE,plot=TRUE){
  ## Iter1: the EM iteration
  ## Iter2: the MH iteration in each EM
  ## burn: the EM iterations for burn-in
  ## MH.burn: the MH interations for burn-in
  ## MH.scale is the step size for MH algorithm
  ## save: whether or not save the intermediate result in the MCEM algorithm
  ## whether or not generate fitting plot 
  Ini.Logratios <- function(base,W){
    ## functions for initialize the logratio responses in LN distribution
    ## base: the base taxa used in logratio transformation
    W.m <- as.matrix(W)
    N <- nrow(W.m)
    Q <- ncol(W.m) #dim of covariates and count vectors
    M=apply(W.m,1,sum)
    
    Psedo_Z.m <- matrix(0,N,Q) 
    Y.m <- matrix(0,N,Q-1) 
    
    for(i in 1:N){
      zeros <- which(W.m[i,]==0)
      nzeros <- which(W.m[i,]!=0)
      Psedo_Z.m[i,zeros] <- (W.m[i,zeros]+0.05)/(M[i]+0.05*length(zeros))
      Psedo_Z.m[i,nzeros]<- (W.m[i,nzeros])/(M[i]+0.05*length(zeros))
      Y.m[i,] <- log(Psedo_Z.m[i,-base]/Psedo_Z.m[i,base])
    }
    attr(Y.m,"center")=apply(Y.m,2,mean)
    attr(Y.m,"base")=paste("the",base,"th taxa")
    Y.m
  }
  
  OLS_fun <- function(X.m,Y.m,b_ini){
    ## function specially for re-estimation
    P <- ncol(X.m); Q <- ncol(Y.m)+1
    CY.m <- scale(Y.m,scale=FALSE)
    b <- matrix(0,P,Q-1)
    for(q in 1:(Q-1)){
      b_q <- b_ini[,q]
      edges <- which(b_q!=0) 
      #the estimated edges for yq
      if(length(edges)>0){
        b_q[edges] <- ginv(t(X.m[,edges])%*%X.m[,edges])%*%t(X.m[,edges])%*%CY.m[,q]
      }
      b[,q]<- b_q
    }
    b
  }
  
  LogisticT <- function(base,b){
    ## to find the perturbation vectors based on the LNM model coefficients
    if(is.vector(b)){
      b.full <- rep(0,length(b)+1)
      b.full[-base] <- b
      b.simplex <- exp(b.full)/sum(exp(b.full))
    }
    
    if(is.matrix(b)){
      b.full <- matrix(0,nrow(b),ncol(b)+1)
      b.full[,-base] <- b
      b.simplex <- t(apply(b.full,1,function(x){exp(x)/sum(exp(x))}))
    }
    b.simplex
  }
  
  X.m <- scale(X);W.m <- as.matrix(W)
  Y.m <- Ini.Logratios(base,W)
  
  N <- nrow(W.m)
  P <- ncol(X.m)
  Q <- ncol(W.m)
  M <- apply(W.m,1,sum)
  
  if(model == "LN"){
    cat("the LN model\n")
    b0=attr(Y.m,"center")
    CY.m <- Y.m-rep(1,N)%*%t(b0)
    b_ini <- remMap(X.m,CY.m,0,lamda)$phi
    master_pred <- which(apply(abs(b_ini),1,sum)!=0)
    
    if(length(master_pred)==0){
      warning("the penalty parameter assigned is too big, tuning parameter procedure suggested")
    }
    
    if(length(master_pred)>0){
      MLE_fit <- lm(CY.m~X.m[,master_pred])
      b <- MLE_fit$coefficients[-1,]
      ## use group lasso maximum likelihood estimator (MLE) hybrid to re-estimate 
      ## check the prediction performance
      eY.m <- as.matrix(X.m[,master_pred])%*% b+rep(1,N)%*%t(b0)
      exp_eY <- exp(eY.m)
      sum_exp_eY <- apply(exp_eY,1,sum) 
      e_Z.m <- matrix(0,N,Q) 
      e_W.m <- matrix(0,N,Q)
      for(i in 1:N){
        e_Z.m[i,-base] <- exp_eY[i,]/(sum_exp_eY[i]+1)
        e_Z.m[i,base] <- 1/(sum_exp_eY[i]+1)
        e_W.m[i,] <- M[i]*e_Z.m[i,]
      }
      MSPE <- mean(apply((e_W.m-W.m)^2,1,sum)/apply(W.m^2,1,sum))
      
      if(plot==TRUE){
        plot(as.vector(sqrt(W.m)),as.vector(sqrt(e_W.m)),
             xlab=expression(sqrt(Observed)),
             ylab=expression(sqrt(Fitted)),main=paste("Fitting Graph"))
      }
      ##Now transfer the estimated constant and coefficients to the constant and coefficients in simplex
      
      b0.simplex <- LogisticT(base,b0)
      names(b0.simplex) <- paste("taxa",1:Q)
      b.simplex <- LogisticT(base,b)
      Magn <- NULL    ##perturbation magnitude
      temp <- diag(1,Q-1)+rep(1,Q-1)%*%t(rep(1,Q-1))
      for(i in 1:length(master_pred)){
        Magn[i] <- sqrt(b[i,]%*%solve(temp)%*%b[i,])
      }
      b.simplex <- cbind(b.simplex,Magn)
      rownames(b.simplex) <- paste("covariate",master_pred)
      colnames(b.simplex) <- c(paste("taxa",1:Q),"Magnitude")
      
      result <- list(lamda,master_pred,b0.simplex,b.simplex,MSPE)
      names(result) <- c("penalty parameter",
                         "Index for master predictors",
                         "Model fitted baseline compositions",
                         "Model fitted perturbation vectors", 
                         "Mean squared prediction error")
      
      print(result)
    }
  }
  
  if(model=="LNM"){
    cat("the LNM model\n")
    cat("Initializing the parameters\n")
    b0 <- attr(Y.m,"center")
    sigma <- var(Y.m)
    U <- chol(solve(sigma))
    CYu.m <- as.matrix(scale((Y.m-rep(1,N)%*%t(b0))%*%U,scale=FALSE)) 
    b_ini <- remMap(X.m,CYu.m,0,lamda)$phi 
    b <- OLS_fun(X.m,Y.m,b_ini)
    master_pred <- which(apply(abs(b),1,sum)!=0)
    eY.m <- X.m%*%b+rep(1,N)%*%t(b0)
    sigma <- var(Y.m-eY.m)
    
    ## store the results in the process
    list.mp <- list(master_pred)
    list.b <- list(b)
    list.b0 <- list(b0)
    list.sigma <- list(sigma)
    list.acceptance <- NULL
    
    ## Start the MCEM algorithm
    for(em in 1:Iter1){
      cat("EM iteration =",em,"\n")
      MH_samples <- list(NULL)
      ## store all the MH samples in one EM step
      Temp.acpt <- 0
      
      for(i in 1:N){
        ## MH independently proceeded for each individual
        ## therefore (N*Iter2) MH iterations in total in each EM step
        YI.store <- matrix(0,Iter2,Q-1)
        YI <- Y.m[i,] 
        WI <- W.m[i,]
        muI <- eY.m[i,]
        
        for(mc in 1:Iter2){
          YI_n <- YI+rnorm(Q-1,0,MH.scale)
          
          Log_pt1 <- sum(WI)*(log(sum(exp(YI))+1)
                              -log(sum(exp(YI_n))+1))
          Log_pt2 <- sum(WI[-base]*(YI_n-YI))
          Log_pt3 <- -0.5*t(YI_n-muI)%*%solve(sigma)%*%(YI_n-muI)
          Log_pt4 <- -0.5*t(YI-muI)%*%solve(sigma)%*%(YI-muI)
          R_ratio <- Log_pt1+Log_pt2+Log_pt3-Log_pt4
          acceptance <- min(1,exp(R_ratio))
          Temp.acpt <- Temp.acpt+acceptance
          ## the acceptance ratio
          temp<-runif(1)
          if(temp < acceptance){YI=YI_n}
          YI.store[mc,]<- YI
        }
        Y.m[i,]<- apply(YI.store[(MH.burn+1):Iter2,],2,mean)
        ## update the Y.m based on the stablized MH samples
        MH_samples <- c(MH_samples,list(YI.store))
      }
      
      list.acceptance <- c(list.acceptance,list(Temp.acpt/(N*Iter2)))
      ## update b0 and sigma with the MH-samples
      SV_sum<-matrix(0,Q-1,Q-1)
      for(mc in (MH.burn+1):Iter2){
        Y.m_mc <- matrix(0,N,Q-1)
        for(i in 1:N){
          Y.m_mc[i,] <- MH_samples[[i+1]][mc,]
        }
        eps_mc <- Y.m_mc-eY.m
        SV_mc <- t(eps_mc)%*%eps_mc
        SV_sum <- SV_sum+SV_mc
      }
      
      sigma <- SV_sum/(N*(Iter2-MH.burn))
      b0 <- apply(Y.m,2,mean)
      
      ##update b
      U <- chol(solve(sigma))
      CY.m <- Y.m-rep(1,N)%*%t(b0)
      CYu.m <- as.matrix(scale(CY.m%*%U,scale=FALSE)) 
      b_ini <- remMap(X.m,CYu.m,0,lamda)$phi
      b <- OLS_fun(X.m,Y.m,b_ini)
      master_pred <- which(apply(abs(b),1,sum)!=0)
      cat("master predictors =",master_pred,"\n")
      
      eY.m <- X.m%*%b+rep(1,N)%*%t(b0)
      list.mp <- c(list.mp,list(master_pred))
      list.b <- c(list.b,list(b))
      list.b0 <- c(list.b0,list(b0))
      list.sigma <- c(list.sigma,list(sigma))    
    }
    
    if(save==TRUE){
      save(list.b,list.b0,list.sigma,list.acceptance,file="listfull.RData")
    }
    
    sigma.sum <- matrix(0,Q-1,Q-1)
    b0.sum <- rep(0,Q-1)
    
    for(em in (burn+1):Iter1){
      sigma.sum <- sigma.sum+list.sigma[[em+1]]
      b0.sum <- b0.sum+list.b0[[em+1]]
    }
    
    sigma.mean <- sigma.sum/(Iter1-burn)
    b0.mean <- b0.sum/(Iter1-burn)
    
    b.med <- matrix(0,P,Q-1)
    # b.mean=matrix(0,P,Q-1)
    for(i in 1:P){
      for(j in 1:(Q-1)){
        bij_store <- NULL
        for(em in (burn+1):Iter1){
          bij_store <- c(bij_store,list.b[[em+1]][i,j])
        }
        b.med[i,j] <- median(bij_store)
        # b.mean[i,j]=mean(b.mean)
      }
    }
    
    master_pred <- which(apply(abs(b.med),1,sum)!=0)
    
    if(length(master_pred)== 0){
      warning("the penalty parameter assigned is too big, tuning parameter procedure suggested")
    }
    
    
    if(length(master_pred)>0){
      eY.m <- X.m%*%b.med+rep(1,N)%*%t(b0.mean)
      exp_eY <- exp(eY.m)
      sum_exp_eY <- apply(exp_eY,1,sum) 
      e_Z.m <- matrix(0,N,Q) 
      e_W.m <- matrix(0,N,Q)
      for(i in 1:N){
        e_Z.m[i,-base] <- exp_eY[i,]/(sum_exp_eY[i]+1)
        e_Z.m[i,base] <- 1/(sum_exp_eY[i]+1)
        e_W.m[i,] <- M[i]*e_Z.m[i,]
      }
      MSPE <- mean(apply((e_W.m-W.m)^2,1,sum)/apply(W.m^2,1,sum))
      
      if(plot==TRUE){
        plot(as.vector(sqrt(W.m)),as.vector(sqrt(e_W.m)),
             xlab=expression(sqrt(Observed)),
             ylab=expression(sqrt(Fitted)),main=paste("Fitting Graph"))
      }
      
      ##Now transfer the estimated constant and coefficients to the constant and coefficients in simplex
      b0.simplex <- LogisticT(base,b0.mean)
      names(b0.simplex) <- paste("taxa",1:Q)
      b.simplex <- LogisticT(base,b.med[master_pred,])
      
      ##perturbation magnitude
      Magn <- NULL
      temp <- diag(1,Q-1)+rep(1,Q-1)%*%t(rep(1,Q-1))
      for(i in 1:length(master_pred)){
        Magn[i] <- sqrt(b.med[master_pred,][i,]%*%solve(temp)%*%b.med[master_pred,][i,])
      }
      b.simplex <- cbind(b.simplex,Magn)
      rownames(b.simplex) <- paste("covariate",master_pred)
      colnames(b.simplex) <- c(paste("taxa",1:Q),"Magnitude")
      
      result <- list(lamda,master_pred,b0.simplex,b.simplex,MSPE)
      names(result) <- c("penalty parameter",
                         "Index for master predictors",
                         "Model fitted baseline compositions",
                         "Model fitted perturbation vectors", 
                         "Mean squared prediction error")
      
      print(result)
    }
  }
}
