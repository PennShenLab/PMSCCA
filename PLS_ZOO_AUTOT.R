PLS_ZOO_AUTOT<-function(X, Y, E,method="PM_SPLS", standardize="Center", c1=NULL,c2=NULL,lambda=0.25,k=5,stratified=NULL,criteria="covariance",seed=0, u0=NULL,v0=NULL, maxiter=100, tolvar=10e-3, tolfun=10e-6, QUIET=T){
  # --------------------------------------------------------------------
  # PM-SPLS algorithm autotune result
  # --------------------------------------------------------------------
  
  # PM-SPLS with Auto-tune (cross validation - training and validation):
  #     This Auto-tune algorithm will automatically tune the algorithm based on the validation canonical correlation;
  #     Input is same as PN_SPLS except we accept vectors of c1, c2 and lambda as hyperparameter candidates;
  #     Each combination of hyperparameter set will be evaluated, one set of hyperparameters which gives the highest
  #       absolute mean of validation canonical correlation will be selected. Then the model will refit the data using
  #       best set of hyperparameters.
  
  # Input:
  #     Required:
  #         X: n-by-p data matrix. Each row of X is a data point. Each column of X is a variable or feature.
  #         Y: n-by-q data matrix. Each row of Y is a data point. Each column of Y is a variable or feature.  
  #         E: p-by-q preference matrix in which all the elements are nonnegative.
  #     Optimal:
  #         k: number of folds for the cross validation (default 4).
  #         stratified: 
  #             NULL - No stratification will be used when separating the dataset into training and validation (default).
  #             length n vector - The folds will be split stratified on the input vector (eg. diagnosis for each patient).
  #         seed: set seed when separating the training and validation dataset (default 0).
  #         c1: vector of regularization parameter (default NULL) # upper bound sqrt(p);
  #         c2: vector of regularization parameter (default NULL) # upper bound sqrt(q);
  #         lambda: vector of regularization parameter (default 0.25).
  #         criteria: 
  #             "covariance" - select model based on the best covariance
  #             "correlation" - select model based on the best correlation
  
  # Author: Jingxuan Bao 
  # Last revision: July 12, 2021
  
  list.of.packages <- c("splitTools")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) {
    install.packages(new.packages)
  }
  library(splitTools)
  if (QUIET==F){
    cat("---------- PM-SPLS AutoTune ----------\n")
    cat("We use",k,"folds cross validation;\n")
    cat("There are",length(lambda),"lambda hyperparameters;\n")
    cat("There are",length(c1),"c1 hyperparameters;\n")
    cat("There are",length(c2),"c2 hyperparameters;\n")
    cat("In total:",(length(lambda)*length(c1)*length(c2)),"hyperparameter settings.\n")
  }
  
  ny <- dim(Y)[1]
  nx <- dim(X)[1]
  if (nx==ny){
    n = nx
  }else{
    stop("Number of subjects for X does not matches number of subjects for Y, please check your input data.")
  }
  p <- dim(X)[2]
  q <- dim(Y)[2]
  
  if (is.null(stratified)){
    idx<-create_folds(1:nrow(X),k = k,type = "basic",shuffle = T,seed = seed)
  }else{
    idx<-create_folds(stratified,k = k,type = "stratified",shuffle = T,seed = seed)
  }
  
  
  
  ## Center and Scale
  if (standardize=="CenterScale"){
    
    X <- X - matrix(apply(X,2,mean), nrow=n, ncol=length(apply(X,2,mean)), byrow=TRUE)
    
    for (i in 1:p) {
      if (!all(X[,i]==0)){
        X[,i] = X[,i]/sd(X[,i])
      }
    }
    
    Y <- Y - matrix(apply(Y,2,mean), nrow=n, ncol=length(apply(Y,2,mean)), byrow=TRUE)
    
    for (j in 1:q) {
      if (!all(Y[,j]==0)){
        Y[,j] = Y[,j]/sd(Y[,j])
      }
    }
    
  }else if(standardize=="Center"){
    X <- X - matrix(apply(X,2,mean), nrow=n, ncol=length(apply(X,2,mean)), byrow=TRUE)
    Y <- Y - matrix(apply(Y,2,mean), nrow=n, ncol=length(apply(Y,2,mean)), byrow=TRUE)
  }else if(standardize=="Raw"){
    warning("The data will not be standardized;\nNote PLS requires all features (i.e., columns) in X, Y to be centered at zero.")
  }else{
    stop("standardize must be one of the following options:\n\"CenterScale\": Center and scale the data;\n(default) \"Center\": Center the data only;\n\"Raw\": No preprocessing will be performed.")
  }
  
  
  
  if (method=="PM_SPLS"){
    if (is.null(c1)){
      warning("Tuning parameters for c1 is not specified, default will be used;")
      c1 <- seq(1,sqrt(ncol(X)),length.out=10)
    }
    if (is.null(c2)){
      warning("Tuning parameters for c2 is not specified, default will be used;")
      c2 <- seq(1,sqrt(ncol(Y)),length.out=10)
    }
    curr_CC<-0
    iter<-0
    for (lambda_num in 1:length(lambda)) {
      for (c1_num in 1:length(c1)) {
        for (c2_num in 1:length(c2)) {
          iter<-iter+1
          if(QUIET==F){
            cat("Hyperparameter Set:",iter,"out of",(length(lambda)*length(c1)*length(c2)),"\n")
          }
          val_res<-rep(0,length(idx))
          for (folds in 1:length(idx)) {
            X_train<-X[idx[[folds]],]
            Y_train<-Y[idx[[folds]],]
            X_train <- X_train - matrix(apply(X_train,2,mean), nrow=nrow(X_train), ncol=length(apply(X_train,2,mean)), byrow=TRUE)
            Y_train <- Y_train - matrix(apply(Y_train,2,mean), nrow=nrow(Y_train), ncol=length(apply(Y_train,2,mean)), byrow=TRUE)
            X_val<-X[!(1:nrow(X))%in%idx[[folds]],] 
            Y_val<-Y[!(1:nrow(Y))%in%idx[[folds]],] 
            X_val <- X_val - matrix(apply(X_train,2,mean), nrow=nrow(X_val), ncol=length(apply(X_train,2,mean)), byrow=TRUE)
            Y_val <- Y_val - matrix(apply(Y_train,2,mean), nrow=nrow(Y_val), ncol=length(apply(Y_train,2,mean)), byrow=TRUE)
            train_res<-PM_SPLS(X_train, Y_train, E, standardize="Raw", c1=c1[c1_num],c2=c2[c2_num],lambda=lambda[lambda_num], u0=u0,v0=v0, maxiter=maxiter, tolvar=tolvar, tolfun=tolfun, QUIET=QUIET,autoT=T)
            # validation data should be centered 
            if(criteria=="correlation"){
              val_res[folds]<-cor(X_val%*%train_res$u,Y_val%*%train_res$v)[1]
            }else if (criteria=="covariance"){
              #val_res[folds]<-t(X_val%*%train_res$u)%*%(Y_val%*%train_res$v)/(norm(train_res$u,type="F") * norm(train_res$v,type="F"))/(nrow(X_val)-1)
              val_res[folds]<-t(X_val%*%train_res$u)%*%(Y_val%*%train_res$v)/( sqrt(sum(train_res$u^2)) * sqrt(sum(train_res$v^2)))/(nrow(X_val)-1)
            }else{
              stop("criteria must be one of \"covariance\" or \"correlation\".")
            }
            
            
          }
          if(abs(mean(val_res, na.rm = T))>curr_CC){
            curr_CC<-abs(mean(val_res))
            cat("Best model is detected: validation canonical correlation =",curr_CC,"\n")
            train_performance<-train_res$cano_corr
            val_performance<-curr_CC # absolute mean of the validation canonical correlation
            best_c1<-c1[c1_num]
            best_c2<-c2[c2_num]
            best_lambda<-lambda[lambda_num]
          }
        }
      }
    }
  }else if(method=="PM_PLS"){
    curr_CC<-0
    iter<-0
    for (lambda_num in 1:length(lambda)) {
      iter<-iter+1
      if(QUIET==F){
        cat("Hyperparameter Set:",iter,"out of",(length(lambda)*length(c1)*length(c2)),"\n")
      }
      val_res<-rep(0,length(idx))
      for (folds in 1:length(idx)) {
        X_train<-X[idx[[folds]],]
        Y_train<-Y[idx[[folds]],]
        X_train <- X_train - matrix(apply(X_train,2,mean), nrow=nrow(X_train), ncol=length(apply(X_train,2,mean)), byrow=TRUE)
        Y_train <- Y_train - matrix(apply(Y_train,2,mean), nrow=nrow(Y_train), ncol=length(apply(Y_train,2,mean)), byrow=TRUE)
        X_val<-X[!(1:nrow(X))%in%idx[[folds]],] 
        Y_val<-Y[!(1:nrow(Y))%in%idx[[folds]],] 
        X_val <- X_val - matrix(apply(X_train,2,mean), nrow=nrow(X_val), ncol=length(apply(X_train,2,mean)), byrow=TRUE)
        Y_val <- Y_val - matrix(apply(Y_train,2,mean), nrow=nrow(Y_val), ncol=length(apply(Y_train,2,mean)), byrow=TRUE)
        train_res<-PM_PLS(X_train, Y_train, E, standardize="Raw",lambda=lambda[lambda_num], u0=u0,v0=v0, maxiter=maxiter, tolvar=tolvar, tolfun=tolfun, QUIET=QUIET,autoT=T)
        # validation data should be centered 
        if(criteria=="correlation"){
          val_res[folds]<-cor(X_val%*%train_res$u,Y_val%*%train_res$v)[1]
        }else if (criteria=="covariance"){
          #val_res[folds]<-t(X_val%*%train_res$u)%*%(Y_val%*%train_res$v)/(norm(train_res$u,type="F") * norm(train_res$v,type="F"))/(nrow(X_val)-1)
          val_res[folds]<-t(X_val%*%train_res$u)%*%(Y_val%*%train_res$v)/( sqrt(sum(train_res$u^2)) * sqrt(sum(train_res$v^2)))/(nrow(X_val)-1)
        }else{
          stop("criteria must be one of \"covariance\" or \"correlation\".")
        }
        
        
      }
      if(abs(mean(val_res, na.rm = T))>curr_CC){
        curr_CC<-abs(mean(val_res))
        cat("Best model is detected: validation canonical correlation =",curr_CC,"\n")
        train_performance<-train_res$cano_corr
        val_performance<-curr_CC # absolute mean of the validation canonical correlation
        best_c1<-NULL
        best_c2<-NULL
        best_lambda<-lambda[lambda_num]
      }
    }
  }else if(method=="SPLS"){
    if (is.null(c1)){
      warning("Tuning parameters for c1 is not specified, default will be used;")
      c1 <- seq(1,sqrt(ncol(X)),length.out=10)
    }
    if (is.null(c2)){
      warning("Tuning parameters for c2 is not specified, default will be used;")
      c2 <- seq(1,sqrt(ncol(Y)),length.out=10)
    }
    curr_CC<-0
    iter<-0
    for (c1_num in 1:length(c1)) {
      for (c2_num in 1:length(c2)) {
        iter<-iter+1
        if(QUIET==F){
          cat("Hyperparameter Set:",iter,"out of",(length(lambda)*length(c1)*length(c2)),"\n")
        }
        val_res<-rep(0,length(idx))
        for (folds in 1:length(idx)) {
          X_train<-X[idx[[folds]],]
          Y_train<-Y[idx[[folds]],]
          X_train <- X_train - matrix(apply(X_train,2,mean), nrow=nrow(X_train), ncol=length(apply(X_train,2,mean)), byrow=TRUE)
          Y_train <- Y_train - matrix(apply(Y_train,2,mean), nrow=nrow(Y_train), ncol=length(apply(Y_train,2,mean)), byrow=TRUE)
          X_val<-X[!(1:nrow(X))%in%idx[[folds]],] 
          Y_val<-Y[!(1:nrow(Y))%in%idx[[folds]],] 
          X_val <- X_val - matrix(apply(X_train,2,mean), nrow=nrow(X_val), ncol=length(apply(X_train,2,mean)), byrow=TRUE)
          Y_val <- Y_val - matrix(apply(Y_train,2,mean), nrow=nrow(Y_val), ncol=length(apply(Y_train,2,mean)), byrow=TRUE)
          train_res<-SPLS(X_train, Y_train, standardize="Raw", c1=c1[c1_num],c2=c2[c2_num], u0=u0,v0=v0, maxiter=maxiter, tolvar=tolvar, tolfun=tolfun, QUIET=QUIET,autoT=T)
          # validation data should be centered 
          if(criteria=="correlation"){
            val_res[folds]<-cor(X_val%*%train_res$u,Y_val%*%train_res$v)[1]
          }else if (criteria=="covariance"){
            #val_res[folds]<-t(X_val%*%train_res$u)%*%(Y_val%*%train_res$v)/(norm(train_res$u,type="F") * norm(train_res$v,type="F"))/(nrow(X_val)-1)
            val_res[folds]<-t(X_val%*%train_res$u)%*%(Y_val%*%train_res$v)/( sqrt(sum(train_res$u^2)) * sqrt(sum(train_res$v^2)))/(nrow(X_val)-1)
          }else{
            stop("criteria must be one of \"covariance\" or \"correlation\".")
          }
          
          
        }
        if(abs(mean(val_res, na.rm = T))>curr_CC){
          curr_CC<-abs(mean(val_res))
          cat("Best model is detected: validation canonical correlation =",curr_CC,"\n")
          train_performance<-train_res$cano_corr
          val_performance<-curr_CC # absolute mean of the validation canonical correlation
          best_c1<-c1[c1_num]
          best_c2<-c2[c2_num]
          best_lambda<-NULL
        }
      }
    }
  }else if(method=="PLS"){
    curr_CC<-0
    val_res<-rep(0,length(idx))
    for (folds in 1:length(idx)) {
      X_train<-X[idx[[folds]],]
      Y_train<-Y[idx[[folds]],]
      X_train <- X_train - matrix(apply(X_train,2,mean), nrow=nrow(X_train), ncol=length(apply(X_train,2,mean)), byrow=TRUE)
      Y_train <- Y_train - matrix(apply(Y_train,2,mean), nrow=nrow(Y_train), ncol=length(apply(Y_train,2,mean)), byrow=TRUE)
      X_val<-X[!(1:nrow(X))%in%idx[[folds]],] 
      Y_val<-Y[!(1:nrow(Y))%in%idx[[folds]],] 
      X_val <- X_val - matrix(apply(X_train,2,mean), nrow=nrow(X_val), ncol=length(apply(X_train,2,mean)), byrow=TRUE)
      Y_val <- Y_val - matrix(apply(Y_train,2,mean), nrow=nrow(Y_val), ncol=length(apply(Y_train,2,mean)), byrow=TRUE)
      train_res<-PLS(X_train, Y_train,standardize="Raw", u0=u0,v0=v0, maxiter=maxiter, tolvar=tolvar, tolfun=tolfun, QUIET=QUIET,autoT=T)
      # validation data should be centered 
      if(criteria=="correlation"){
        val_res[folds]<-cor(X_val%*%train_res$u,Y_val%*%train_res$v)[1]
      }else if (criteria=="covariance"){
        #val_res[folds]<-t(X_val%*%train_res$u)%*%(Y_val%*%train_res$v)/(norm(train_res$u,type="F") * norm(train_res$v,type="F"))/(nrow(X_val)-1)
        val_res[folds]<-t(X_val%*%train_res$u)%*%(Y_val%*%train_res$v)/( sqrt(sum(train_res$u^2)) * sqrt(sum(train_res$v^2)))/(nrow(X_val)-1)
      }else{
        stop("criteria must be one of \"covariance\" or \"correlation\".")
      }
      
      
    }
    if(abs(mean(val_res, na.rm = T))>curr_CC){
      curr_CC<-abs(mean(val_res))
      cat("Best model is detected: validation canonical correlation =",curr_CC,"\n")
      train_performance<-train_res$cano_corr
      val_performance<-curr_CC # absolute mean of the validation canonical correlation
      best_c1<-NULL
      best_c2<-NULL
      best_lambda<-NULL
    }
  }else{
    stop("method should be one of PM_SPLS, PM_PLS, SPLS, PLS!")
  }
  
  
  if (method=="PM_SPLS"){
    refit_res<-PM_SPLS(X, Y, E, standardize="Raw", c1=best_c1,c2=best_c2,lambda=best_lambda, u0=u0,v0=v0, maxiter=maxiter, tolvar=tolvar, tolfun=tolfun, QUIET=QUIET,autoT=T)
  }else if(method=="PM_PLS"){
    refit_res<-PM_PLS(X, Y, E, standardize="Raw", lambda=best_lambda, u0=u0,v0=v0, maxiter=maxiter, tolvar=tolvar, tolfun=tolfun, QUIET=QUIET,autoT=T)
  }else if(method=="SPLS"){
    refit_res<-SPLS(X, Y, standardize="Raw", c1=best_c1,c2=best_c2, u0=u0,v0=v0, maxiter=maxiter, tolvar=tolvar, tolfun=tolfun, QUIET=QUIET,autoT=T)
  }else if(method=="PLS"){
    refit_res<-PLS(X, Y, standardize="Raw", u0=u0,v0=v0, maxiter=maxiter, tolvar=tolvar, tolfun=tolfun, QUIET=QUIET,autoT=T)
  }else{
    stop("method should be one of PM_SPLS, PM_PLS, SPLS, PLS!")
  }
  
  return(list(bestModel=refit_res,training_performance=train_performance,validation_performance=val_performance,best_lambda=best_lambda,best_c1=best_c1,best_c2=best_c2))
}


