PLS<-function(X, Y,  standardize="Center", u0=NULL,v0=NULL, maxiter=100, tolvar=10e-3, tolfun=10e-6, QUIET=F,autoT=F){
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
  
  # Author: Jingxuan Bao, Jiahang Sha, Kefei Liu 
  # Last revision: Aug 1, 2022
  
  
  output<-PM_PLS(X, Y, E=matrix(0,nrow = ncol(X),ncol = ncol(Y)), standardize=standardize, lambda=0, u0=u0,v0=v0, maxiter=maxiter, tolvar=tolvar, tolfun=tolfun, QUIET=QUIET,autoT=autoT)
  
  return(output)
}




















