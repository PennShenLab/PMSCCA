PM_SPLS<-function(X, Y, E, standardize="Center", c1=NULL,c2=NULL,lambda=0.25, u0=NULL,v0=NULL, maxiter=100, tolvar=10e-3, tolfun=10e-6, QUIET=F,autoT=F){
  # --------------------------------------------------------------------
  # PM-SPLS algorithm 
  # --------------------------------------------------------------------
  
  # PM-SPLS model:
  #     maximize     (1-lambda) * 1/(n-1) * u' * X' * Y * v + lambda * |u|' * E * |v| 
  #     subject to   ||u||_2^2 <= 1, ||u||_1 <= c_1
  #                  ||v||_2^2 <= 1, ||v||_1 <= c_2
  #     where n is the sample size.
  
  # Input:
  #     Required:
  #         X: n-by-p data matrix. Each row of X is a data point. Each column of X is a variable or feature.
  #         Y: n-by-q data matrix. Each row of Y is a data point. Each column of Y is a variable or feature.  
  #         E: p-by-q preference matrix in which all the elements are nonnegative.
  #     Optional:
  #         standardize: One of the following options 
  #             - "Center" (default): Each feature (i.e., column) in X, Y will be centered at zero;
  #             - "centerScale": Each feature (i.e., column) in X, Y will be standardized to have a mean of zero and a standard deviation of one;
  #             - "Raw": The features (i.e., columns) in X, Y will not be standardized (neither centered nor scaled).
  #         Note that according to the definition, PLS requires all features (i.e., columns) in X, Y to be centered at zero. 
  #         Thus, unless the data have already been centered before calling the PLS function, one generally should not set center = false. 
  
  #         c1: regularization parameter, range [0,1);
  #         c2: regularization parameter, range [0,1);
  #         lambda: regularization parameter, range [0,1).
  
  #         u0: initial guess for u. 
  #         v0: initial guess for v. 
  
  #         maxiter: maximum number of iterations allowed. The default value is 100.
  #         tolvar: tolerance level for the variables. The default value is 10^{-3}.
  #         tolfuc: tolerance level for the objective function value. The default value is 10^{-6}.
  #         QUIET: (true/false): used to monitor the running of the algorithm for diagnostics, reporting, termination checks. 
  #                If QUIET = false (default), intermediate results at each iteration will appear on the screen.
  
  # Output:
  #     u: p-by-1 weight vector for data X.
  #     v: q-by-1 weight vector for data Y.
  #     objval: objective function value.
  #     cano_cov: canonical covariance.
  #     cano_corr: canonical correlation.
  #     d: a constant used to calculate the residual covariance matrices and subsequent canonical components.
  #     history: a structure that contains the objective and constraint values at each iteration and the differences between successive iterations.

  # Note:
  #     This code is adapted from Kefei's MatLab code
  
  # References:
  # 1. Kefei Liu, Long Qi and Li Shen, "Enhanced canonical correlation analysis model via elastic net regularization".
  # 2. Kefei Liu and Li Shen, "Some notes on sparse canonical correlation analysis".
  
  # Author: Jingxuan Bao, Kefei Liu 
  # Last revision: July 4, 2021
  

  ## Package checking
  list.of.packages <- c("pracma")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) {
    install.packages(new.packages)
  }
  library(pracma)
  
  ## Initialize the list to store the means and sds for X and Y before and after preprocess
  parasStandardize <- vector(mode = "list",3)
  # Store the statistics for the oriiginal data
  names(parasStandardize) <- c("before.preprocess","after.preprocess","init.params")
  parasStandardize[["before.preprocess"]] <- list(mean.X = apply(X,2,mean),
                                                  sd.X = apply(X,2,sd),
                                                  mean.Y = apply(Y,2,mean),
                                                  sd.Y = apply(Y,2,sd))
  
  
  ## Input check
  ny <- dim(Y)[1]
  nx <- dim(X)[1]
  if (nx==ny){
    n = nx
    cat("Number of subjects is",n,"\n")
  }else{
    stop("Number of subjects for X does not matches number of subjects for Y, please check your input data.")
  }
  p <- dim(X)[2]
  q <- dim(Y)[2]
  cat("Number of features for X is",p,"\n")
  cat("Number of features for Y is",q,"\n")
  
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
  
  # only center
  }else if(standardize=="Center"){
    X <- X - matrix(apply(X,2,mean), nrow=n, ncol=length(apply(X,2,mean)), byrow=TRUE)
    Y <- Y - matrix(apply(Y,2,mean), nrow=n, ncol=length(apply(Y,2,mean)), byrow=TRUE)
  
  # no center no scale
  }else if(standardize=="Raw"){
    if(autoT==F){
      warning("The data will not be standardized;\nNote PLS requires all features (i.e., columns) in X, Y to be centered at zero.")
    }
    
  }else{
    stop("standardize must be one of the following options:\n\"CenterScale\": Center and scale the data;\n(default) \"Center\": Center the data only;\n\"Raw\": No preprocessing will be performed.")
  }
  # Store the statistics for the preprocessed data
  parasStandardize[["after.preprocess"]] <- list(mean.X = apply(X,2,mean),
                                                 sd.X = apply(X,2,sd),
                                                 mean.Y = apply(Y,2,mean),
                                                 sd.Y = apply(Y,2,sd))

  
  ## c1, c2, lambda
  if (is.null(c1)){
    cat("The regularization parameter c1 is not provided! (default will be used)\n");
    c1_min <- 1.0
    c1_max <- sqrt(p)
    c1 <- sqrt(c1_min*c1_max)
  }
  if (is.null(c2)){
    cat("The regularization parameter c2 is not provided! (default will be used)\n");
    c2_min <- 1.0
    c2_max <- sqrt(q)
    c2 <- sqrt(c2_min*c2_max)
  }
  
  ## u0, v0
  if (is.null(u0)&&is.null(v0)){
    svd_val<-svds((t(X)%*%Y)/(n-1),1)
    u <- svd_val$u
    v <- svd_val$v
  }else if(is.null(u0)&&!is.null(v0)){
    v <- v0
    u <- t(X)%*%(Y%*%v0)/(n-1);
    u <- u / sqrt(sum(u^2));
  }else if(!is.null(u0)&&is.null(v0)){
    u <- u0
    v <- t(Y)%*%(X%*%u0)/(n-1);
    v <- v / sqrt(sum(v^2));
  }else{
    u <- u0
    v <- v0
  }
  
  if(QUIET==F){
    cat("Parameters: c1 =",c1,"c2 =",c2,"lambda =",lambda,"\n")
  }
  
  ## Store the initial parameters
  parasStandardize[["init.params"]] <- list(c1 = c1,
                                            c2 = c2,
                                            lambda = lambda,
                                            u0 = u,
                                            v0 = v)
  
  
  ## Store the algorithm
  historyMat <- data.frame(matrix(NA,nrow = maxiter+1,ncol = 9))
  colnames(historyMat) <- c('time','L2norm_u','L1norm_u','L2norm_v','L1norm_v','objval','diff_u','diff_v','diff_objval');
  # Reporting
  historyMat[1,"time"] <- format(Sys.time(), "%b %d %X %Y")
  # Debugging and diagnostics
  historyMat[1,"objval"] <- (1-lambda) * (t(X%*%u)%*%(Y%*%v)/(n-1)) + lambda * (abs(t(u)) %*% E %*% abs(v));
  historyMat[1,"L2norm_u"] = sum(u^2, na.rm = T);
  historyMat[1,"L2norm_v"] = sum(v^2, na.rm = T);
  historyMat[1,"L1norm_u"] = sum(abs(u), na.rm = T);
  historyMat[1,"L1norm_v"] = sum(abs(v), na.rm = T);
  # Termination Checkes
  historyMat[1,"diff_u"] = 1;
  historyMat[1,"diff_v"] = 1;
  historyMat[1,"diff_objval"] = historyMat[1,"objval"]
  
  if(QUIET==F){
    cat("START Time:",historyMat[1,"time"],
        "\nL2norm_u =",historyMat[1,"L2norm_u"],
        "L2norm_v =",historyMat[1,"L2norm_v"],
        "\nL1norm_u =",historyMat[1,"L1norm_u"],
        "L1norm_v =",historyMat[1,"L1norm_v"],
        "\n")
  }
  
  ## Algorithm Procedure
  iterNum <- 1
  # while (iterNum<=maxiter&&
  #        (max(historyMat[iterNum,"diff_u"],historyMat[iterNum,"diff_v"])>tolvar||
  #         historyMat[iterNum,"diff_objval"]>tolfun)) {
    
  while (iterNum<=maxiter&&
         (historyMat[iterNum,"diff_u"]+historyMat[iterNum,"diff_v"]>=2*tolvar^2 || 
          abs(historyMat[iterNum,"diff_objval"])>= tolfun)) {
    
    
    # Update u 
    u_old = u;
    a = (1-lambda) * t(X)%*%(Y%*%v)/(n-1);
    b = lambda * E %*% abs(v);
    argu = a + b * sign_binary(a);
    u = QCLP_solver(argu,c1);
    # Done updating u 
    
    
    # Update v %
    v_old = v;
    alpha = (1-lambda) * t(Y)%*%(X%*%u)/(n-1);
    beta = lambda * t(E) %*% abs(u);
    argv = alpha + beta*sign_binary(alpha);
    v = QCLP_solver(argv,c2);
    # Done updating v 
    
    ## Reporting
    historyMat[iterNum+1,"time"] <- format(Sys.time(), "%b %d %X %Y")
    # Debugging and diagnostics
    historyMat[iterNum+1,"objval"] <- (1-lambda) * (t(X%*%u)%*%(Y%*%v)/(n-1)) + lambda * (abs(t(u)) %*% E %*% abs(v));
    historyMat[iterNum+1,"L2norm_u"] = sum(u^2, na.rm = T);
    historyMat[iterNum+1,"L2norm_v"] = sum(v^2, na.rm = T);
    historyMat[iterNum+1,"L1norm_u"] = sum(abs(u), na.rm = T);
    historyMat[iterNum+1,"L1norm_v"] = sum(abs(v), na.rm = T);
    # Termination Checkes
    historyMat[iterNum+1,"diff_u"] = min(sum((u-u_old)^2),sum((u+u_old)^2));
    historyMat[iterNum+1,"diff_v"] = min(sum((v-v_old)^2),sum((v+v_old)^2));
    historyMat[iterNum+1,"diff_objval"] = historyMat[iterNum+1,"objval"] - historyMat[iterNum,"objval"]
    
    
    if(QUIET==F){
      cat("Iteration:",iterNum,"out of",maxiter,"\n")
      # cat("Time:",historyMat[iterNum+1,"time"],
      #     "\nL2norm_u =",historyMat[iterNum+1,"L2norm_u"],
      #     "L2norm_v =",historyMat[iterNum+1,"L2norm_v"],
      #     "\nL1norm_u =",historyMat[iterNum+1,"L1norm_u"],
      #     "L1norm_v =",historyMat[iterNum+1,"L1norm_v"],
      #     "\n")
      # cat("Current objective function value:",historyMat[iterNum+1,"objval"],"\n")
      # cat("Difference of u:",historyMat[iterNum+1,"diff_u"],
      #     "Difference of v:",historyMat[iterNum+1,"diff_v"],
      #     "\nDifference of objective function value:",historyMat[iterNum+1,"diff_objval"],
      #     "\n")
    }
    
    iterNum = iterNum+1;
  }
  
  # Update f 
  objval <- (1-lambda) * (t(X%*%u)%*%(Y%*%v)/(n-1)) + lambda * (abs(t(u)) %*% E %*% abs(v));
  if (any(u!=0) && any(v!=0)){
    cano_cov = t(X%*%u)%*%(Y%*%v)/(n-1) / sqrt(sum(u^2)) / sqrt(sum(v^2));
    cano_corr = t(X%*%u)%*%(Y%*%v) / sqrt(sum((X%*%u)^2)) / sqrt(sum((Y%*%v)^2));
    d = t(X%*%u)%*%(Y%*%v)/(n-1) / sum(u^2) / sum(v^2);
  }else{
    cano_cov = 0;
    cano_corr = 0;
    d = 0;
  }
  # Done updating f
  
  if(QUIET==F){
    if (max(historyMat[iterNum,"diff_u"],historyMat[iterNum,"diff_v"])<=tolvar||
        historyMat[iterNum,"diff_objval"]<=tolfun){
      cat("Algorithm Converges\n")
    }else{
      cat("Algorithm does not converge, maximum iteration is reached.\n")
    }
    
    cat("END Time:",historyMat[iterNum,"time"],
        "\nCanonical Covariance =",cano_cov,
        "\nCanonical Correlation =",cano_corr,
        "\nd Constant =",d,
        "\n")
  }
  
  return(list(u=u, v=v, objval=objval, cano_cov=cano_cov, cano_corr=cano_corr, d=d, history=historyMat[!apply(is.na(historyMat),1,any),], parasStandardize=parasStandardize))
}

sign_binary<-function(u){
  f <- (u>0)-(u<0);
  idx <- which(f==0,arr.ind = T)
  if(length(idx)!=0){
    f[which(f==0)] <- 2*sample(c(0,1),length(which(f==0)),replace = T)-1
  }
  return(f)
}

norm_matlab<-function(A){
  if(1%in%dim(A)){
    normA<-norm(A,"F")
  }else{
    normA<-max(svd(A)$d)
  }
  return(normA)
}

QCLP_solver<-function(a,c){
  p = length(a)
  if (all(a==0)){
    u0 = matrix(1,p,1);
    scale_const = max(sqrt(p),p/c);
    uhat = u0 / scale_const;
  }else if(c >= sum(abs(a/norm_matlab(a)))){
    uhat = a/norm_matlab(a);
  }else{
    S = which(abs(a) == max(abs(a)));
    if (c <= sqrt(length(S))){
      uhat = matrix(0,p,1);
      uhat[S] = c/length(S) * sign(a[S]);
    }else{
      delta = binarySearch(a,c);
      su = shrinkage(a,delta);
      uhat = su/norm_matlab(su);
    }
  }
  return(uhat)
}

myFunc<-function(x,z,c){
  return(sum(z - x)/sqrt(sum((z-x)^2)) - c)
}

binarySearch<-function(a,c){
  a_abs_unique_sored = sort(unique(rbind(0,abs(a))));
  idx1 = 1; 
  idx2 = length(a_abs_unique_sored);
  while((idx2-idx1)>=2){
    idx = floor((idx1+idx2)/2);
    su = shrinkage(a,a_abs_unique_sored[idx]);
    if (sum(abs(su/norm_matlab(su)))<c){
      idx2 = idx;
    }else if (sum(abs(su/norm_matlab(su)))>c){
      idx1 = idx;
    }else{
      delta = a_abs_unique_sored[idx];
      return(delta)
    }
  }
  lam1 = a_abs_unique_sored[idx1];
  lam2 = a_abs_unique_sored[idx2];
  a_abs = abs(a);
  a_abs_thresholded = a_abs[a_abs > lam1];
  delta = fzero(myFunc, c(lam1,lam2), maxiter = 100, tol = .Machine$double.eps^(1/2), z=a_abs_thresholded,c=c)$x;
  return(delta)
}

shrinkage<-function(x, delta){
  return(apply(x - delta,1,max,0) - apply(-x - delta,1,max,0))
}





















