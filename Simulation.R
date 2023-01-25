  # Main file to execute for reproducing simulation studies. 
  # Please change the file directory as appropriate.
  # Author: Bao & Sha
  

rm(list=ls())
source("./PM_SPLS.R")
source("./SPLS.R")
source("./PM_PLS.R")
source("./PLS.R")
source("./PLS_ZOO_AUTOT.R")
library(pracma)
library(rARPACK)

## generate E
sub2ind<-function(sizeA,r,c){
  nr = sizeA[1]
  nc = sizeA[2]
  return((c-1)*nr + r)
}
generate_preference<-function(gene_NUM, img_NUM, gene_IDX_relevant, img_IDX_relevant,mean_logfc=1, sd_logfc=0.5,pcnt_relevant_diffExpr=0.6, ratio_nonrelevant_occupy=0.05){
  list.of.packages <- c("truncnorm","pracma")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) {
    install.packages(new.packages)
  }
  library(truncnorm)
  
  E = matrix(0,nrow=gene_NUM,ncol=img_NUM)
  
  # Relevant block expression pattern
  gene_NUM_relevant = length(gene_IDX_relevant)
  img_NUM_relevant = length(img_IDX_relevant)
  
  E_relevant = matrix(0,nrow=gene_NUM_relevant,ncol=img_NUM_relevant)
  NUM_relevant_w_diffExpr = round(gene_NUM_relevant * img_NUM_relevant * pcnt_relevant_diffExpr)
  E_relevant[sample(1:(gene_NUM_relevant*img_NUM_relevant), size=NUM_relevant_w_diffExpr)] = rtruncnorm(NUM_relevant_w_diffExpr,a=0,b=Inf,mean =mean_logfc,sd=sd_logfc)
  
  
  E[gene_IDX_relevant, img_IDX_relevant] =  E_relevant
  
  
  # Non-relevant blocks of expression pattern
  n_nonrelevant = floor((gene_NUM*img_NUM - gene_NUM_relevant*img_NUM_relevant) * ratio_nonrelevant_occupy)
  
  temp = meshgrid(img_IDX_relevant,gene_IDX_relevant)
  img_Sub_relevant = temp$X
  gene_Sub_relevant = temp$Y
  
  IND_nonrelevant = sample(setdiff(1:(gene_NUM*img_NUM), sub2ind(c(gene_NUM,img_NUM), as.vector(gene_Sub_relevant), as.vector(img_Sub_relevant))), size=n_nonrelevant)
  
  E[IND_nonrelevant] = rtruncnorm(n_nonrelevant,a=0,b=Inf,mean = mean_logfc,sd=sd_logfc) 
  
  return(E)
}



## Paramter set up
n = 2500  # number of subjects
p = 200 # number of X variables
q = 100  # number of Y variables

X = as.matrix(read.csv("C:\\Users\\shenl\\Box\\cca notes\\xrho031.txt", header = FALSE))
Y = as.matrix(read.csv("C:\\Users\\shenl\\Box\\cca notes\\yrho031.txt", header = FALSE))
C = as.matrix(read.csv("C:\\Users\\shenl\\Box\\cca notes\\crho031.txt", header = FALSE))
Xvar_weight = as.matrix(read.csv("C:\\Users\\shenl\\Box\\cca notes\\vrho031.txt", header = FALSE))
Yvar_weight = as.matrix(read.csv("C:\\Users\\shenl\\Box\\cca notes\\wrho031.txt", header = FALSE))


# # relevant genetic features
# NUM_Xvar_relevant = 20                  # number of relevant genetic features
# IDX_Xvar_relevant = 1:NUM_Xvar_relevant # location of relevant genetic features
# # IDX_Xvar_relevant = sort(datasample(1:p,NUM_Xvar_relevant,'Replace',false))
# 
# # relevant imaging features
# NUM_Yvar_relevant = 30                  # number of relevant imaging features
# IDX_Yvar_relevant = 1:NUM_Yvar_relevant  # location of relevant imaging features
# # IDX_Yvar_relevant = sort(datasample(1:q,NUM_Yvar_relevant,'Replace',false)) # location of relevant imaging features


#rho = 0.8 #0.375#0.8#0.25 # target canonical correlation coefficient !!!NOT ACTUALLY USED!!!
IDX_Xvar_relevant_A = 31:50
IDX_Xvar_relevant_B = 61:90
IDX_Yvar_relevant_A = 21:30
IDX_Yvar_relevant_B = 71:90
# 
# # Generate latent variable that captures the variation shared by both data sets
# z = rnorm(n)
# 
# ## Generate data X
# Xvar_weight = rep(0,p)
# Xvar_weight[IDX_Xvar_relevant] = 1
# sigma_X_sq = (1-rho)/rho * sum(Xvar_weight^2)
# X = z %*% t(Xvar_weight) + sqrt(sigma_X_sq) * matrix(rnorm(n*p),nrow = n,ncol = p)


# ## Generate data Y
# Yvar_weight = rep(0,q)
# Yvar_weight[IDX_Yvar_relevant] = 1
# sigma_Y_sq = (1-rho)/rho * sum(Yvar_weight^2)
# Y = z %*% t(Yvar_weight) + sqrt(sigma_Y_sq) * matrix(rnorm(n*p),nrow = n,ncol = q)
# 
# 
# ## Population covariance matrices
# Sigma_XY = Xvar_weight %*% t(Yvar_weight)
# Sigma_XX = Xvar_weight %*% t(Xvar_weight) + sigma_X_sq * diag(p)
# Sigma_YY = Yvar_weight %*% t(Yvar_weight) + sigma_Y_sq * diag(q)
# 
# ## CCA theoretical parameters
# ## Canonical weights
# CCA.theo.u = Xvar_weight / sqrt(sum(Xvar_weight^2)) / sqrt(sum(Xvar_weight^2)+sigma_X_sq)
# CCA.theo.v = Yvar_weight / sqrt(sum(Yvar_weight^2)) / sqrt(sum(Yvar_weight^2)+sigma_Y_sq)
# 
# CCA.theo.cano_cov = sqrt( sum(Xvar_weight^2) * sum(Yvar_weight^2) )
# CCA.theo.cano_corr = sqrt( sum(Xvar_weight^2) / (sum(Xvar_weight^2)+sigma_X_sq) * sum(Yvar_weight^2) / (sum(Yvar_weight^2)+sigma_Y_sq) )
# 
# ## PLS theoretical parameters
# ## Canonical weights/loadings
PLS.theo.u = Xvar_weight / sqrt(sum(Xvar_weight^2))
PLS.theo.v = Yvar_weight / sqrt(sum(Yvar_weight^2))

# PLS.theo.cano_cov = sqrt( sum(Xvar_weight^2) * sum(Yvar_weight^2) )
# PLS.theo.cano_corr = sqrt( sum(Xvar_weight^2) / (sum(Xvar_weight^2)+sigma_X_sq) * sum(Yvar_weight^2) / (sum(Yvar_weight^2)+sigma_Y_sq) )

E = generate_preference( p, q, IDX_Xvar_relevant_A, IDX_Yvar_relevant_A, 1, 0.5, 0.4, 0.02 )+
  generate_preference( p, q, IDX_Xvar_relevant_B, IDX_Yvar_relevant_A, 1, 0.5, 0.3, 0.0 )+
  generate_preference( p, q, IDX_Xvar_relevant_A, IDX_Yvar_relevant_B, 1, 0.5, 0.3, 0.01 )+
  generate_preference( p, q, IDX_Xvar_relevant_B, IDX_Yvar_relevant_B, 1, 0.5, 0.4, 0.02 )
plot_ly(z = E , type = "heatmap")

X_centered = X - matrix(apply(X,2,mean), nrow=n, ncol=length(apply(X,2,mean)), byrow=TRUE)
Y_centered = Y - matrix(apply(Y,2,mean), nrow=n, ncol=length(apply(Y,2,mean)), byrow=TRUE)
# E = E / norm(E,"F") * norm(t(X_centered)%*%Y_centered / (nrow(X_centered)-1),"F")


## Set hyperparameters
# lambda_set = c(0, 0.0001,0.0005, 0.001,0.005, 0.01,0.05, 0.1, 0.25, 0.5)
# lambda_set = sort(unique(c(lambda_set,1-lambda_set)))
lambda_set = seq(0.4,0.95,length =20)
c1 <- seq(1,sqrt(ncol(X)),length.out=20)

c2 <- seq(1,sqrt(ncol(Y)),length.out=20)


## Start simulation
## Plot
# Libraries
library(ggplot2)
plot_lollipop<-function(train_model,method,weight){
  methods_name<-c("PM_PLS","SPLS","PLS","PM_SPLS")
  if(!method%in%methods_name){
    stop("method should be one of PM_SPLS, PM_PLS, SPLS, PLS!")
  }
  auto_t_path<-file.path("./simulation", method)
  if (!dir.exists(auto_t_path)){
    dir.create(auto_t_path)
  }
  if (weight=="u"){
    ## u
    # Create data
    data.u <- data.frame(
      feature=c(1:ncol(X),1:ncol(X)),
      u=c(PLS.theo.u,train_model$bestModel$u),
      estimated=c(rep(0,ncol(X)),rep(1,ncol(X)))
    )
    
    theme_new <- theme_set(theme_minimal(base_size = 40))
    theme_new <- theme_update(legend.position="bottom")
    lollipop<-ggplot(data.u, aes(y = feature, x = u, label = u, fill = estimated, colour = estimated)) +
      geom_segment(aes(x = 0, y = feature, xend = u, yend = feature), color = "grey50", size = 0.75) +
      geom_point(size = 3) +
      facet_wrap(~estimated) +
      ylab("Feature") + 
      xlab("Weights") + 
      ggtitle("Lollipop Plot for X Weight (True vs Estimated)") + 
      guides(fill = "none",color = "none",scale="none")
    
    png(filename = paste0("./simulation/",method,"/Simulation_u.png"),width = 1400,height = 1400,bg ="white")
    print(lollipop)
    dev.off()
  }else if (weight=="v"){
    ## v
    # Create data
    data.v <- data.frame(
      feature=c(1:ncol(Y),1:ncol(Y)),
      v=c(PLS.theo.v,train_model$bestModel$v),
      estimated=c(rep(0,ncol(Y)),rep(1,ncol(Y)))
    )
    
    theme_new <- theme_set(theme_minimal(base_size = 40))
    theme_new <- theme_update(legend.position="bottom")
    lollipop<-ggplot(data.v, aes(y = feature, x = v, label = v, fill = estimated, colour = estimated)) +
      geom_segment(aes(x = 0, y = feature, xend = v, yend = feature), color = "grey50", size = 0.75) +
      geom_point(size = 3) +
      facet_wrap(~estimated) +
      ylab("Feature") + 
      xlab("Weights") + 
      ggtitle("Lollipop Plot for Y Weight (True vs Estimated)") + 
      guides(fill = "none",color = "none",scale="none")
    
    png(filename = paste0("./simulation/",method,"/Simulation_v.png"),width = 1400,height = 1400,bg ="white")
    print(lollipop)
    dev.off()
  }else{
    stop("weight must be u or v!")
  }
  return(lollipop)
}

library(splitTools)
split_idx<-create_folds(1:nrow(X),k = 5,type = "basic",shuffle = T,seed = 20)
X_train_val<-X[split_idx[[1]],]
Y_train_val<-Y[split_idx[[1]],]
X_test<-X[!(1:nrow(X))%in%split_idx[[1]],] 
Y_test<-Y[!(1:nrow(Y))%in%split_idx[[1]],] 

X_test <- X_test - matrix(apply(X_train_val,2,mean), nrow=nrow(X_test), ncol=length(apply(X_train_val,2,mean)), byrow=TRUE)
Y_test <- Y_test - matrix(apply(Y_train_val,2,mean), nrow=nrow(Y_test), ncol=length(apply(Y_train_val,2,mean)), byrow=TRUE)
E = E / svds(E,1)$d * svds(t(X_train_val)%*%Y_train_val,1)$d / (nrow(X_train_val)-1)
library(plotly)
plot_ly(z=E, type="heatmap")
res_table<-data.frame(matrix(NA,nrow = 1,ncol = 4))
# methods_name<-c("PM_PLS","PLS","SPLS","PM_SPLS")
methods_name<-c("PM_PLS","PLS")
colnames(res_table)<-methods_name
for (idx_method in 1:length(methods_name)) {
  auto_t_path<-file.path("./simulation", methods_name[idx_method])
  if (!dir.exists(auto_t_path)){
    dir.create(auto_t_path)
  }
  if (methods_name[idx_method]=="PM_SPLS"){
    train_model <- PLS_ZOO_AUTOT(X_train_val, 
                                 Y_train_val, 
                                 E,method = methods_name[idx_method], 
                                 c1=c1, c2=c2,
                                 lambda = lambda_set,
                                 # c1=c(5.84),c2=c(6.68),
                                 # lambda=c(5e-4),
                                 k=5)
    # plot_lollipop(train_model,"PM_SPLS","u")
    # plot_lollipop(train_model,"PM_SPLS","v")
  }else if(methods_name[idx_method]=="PM_PLS"){
    train_model <- PLS_ZOO_AUTOT(X_train_val, 
                                 Y_train_val, 
                                 E,method = methods_name[idx_method], 
                                 # lambda=c(0.05),
                                 lambda=lambda_set,
                                 k=5)
    # plot_lollipop(train_model,"PM_PLS","u")
    # plot_lollipop(train_model,"PM_PLS","v")
  }else if(methods_name[idx_method]=="SPLS"){
    train_model <- PLS_ZOO_AUTOT(X_train_val, 
                                 Y_train_val, 
                                 E,method = methods_name[idx_method], 
                                 # c1=c(5.8418),c2=c(1.4736),
                                 c1=c1,c2=c2,
                                 k=5)
    # plot_lollipop(train_model,"SPLS","u")
    # plot_lollipop(train_model,"SPLS","v")
  }else if(methods_name[idx_method]=="PLS"){
    train_model <- PLS_ZOO_AUTOT(X_train_val, 
                                 Y_train_val, 
                                 E,method = methods_name[idx_method], 
                                 k=5)
    # plot_lollipop(train_model,"PLS","u")
    # plot_lollipop(train_model,"PLS","v")
  }
  
  saveRDS(train_model,file = paste0(auto_t_path,"/best_train0.rds"))
  saveRDS(X,file = paste0(auto_t_path,"/Simu_X.rds"))
  saveRDS(Y,file = paste0(auto_t_path,"/Simu_Y.rds"))
  saveRDS(E,file = paste0(auto_t_path,"/Simu_E.rds"))
  saveRDS(X_train_val,file = paste0(auto_t_path,"/Simu_X_train_val.rds"))
  saveRDS(Y_train_val,file = paste0(auto_t_path,"/Simu_Y_train_val.rds"))
  saveRDS(X_test,file = paste0(auto_t_path,"/Simu_X_test.rds"))
  saveRDS(Y_test,file = paste0(auto_t_path,"/Simu_Y_test.rds"))
  ## Testing performance
  res_table[1,idx_method]<-cor(X_test%*%train_model$bestModel$u,Y_test%*%train_model$bestModel$v)[1]
}
 write.table(res_table,file = "./simulation/testing_corr.csv",sep = ",",col.names = T,row.names = F,quote = F)



BT = readRDS("C:\\Users\\shenl\\OneDrive\\Desktop\\0814 cov + vali jiaozheng\\simulation\\PM_PLS\\best_train0.rds")
BTpls = readRDS("C:\\Users\\shenl\\OneDrive\\Desktop\\0814 cov + vali jiaozheng\\simulation\\PLS\\best_train0.rds")
pmu = as.matrix(BT$bestModel$u)
pmv = as.matrix(BT$bestModel$v)
plot_ly(z= t(X_train_val) %*% Y_train_val, type="heatmap")
plot_ly(z= pmu %*% t(pmv), type="heatmap")
plsu = as.matrix(BTpls$bestModel$u)
plsv = as.matrix(BTpls$bestModel$v)
temp  = eye(p) - (pmu%*%t(pmu)%*%t(X_train_val)%*% X_train_val / norm(X_train_val%*%pmu, type="2")^2)
Xonepm = X_train_val %*% temp
Yonepm = Y_train_val %*% (pmv%*%t(pmv)%*%t(Y_train_val)%*% Y_train_val / norm(Y_train_val%*%pmv, type="2")^2)
Yonepm = Y_train_val - Yonepm
plot_ly(z= t(Xonepm) %*% Yonepm, type="heatmap")
Xonepls = X_train_val %*% (diag(p) - plsu%*%t(plsu)%*%t(X_train_val)%*% X_train_val / norm(X_train_val%*%plsu, type="2")^2)
Yonepls = Y_train_val %*% (diag(q) - plsv%*%t(plsv)%*%t(Y_train_val)%*% Y_train_val / norm(Y_train_val%*%plsv, type="2")^2)
plot_ly(z= t(Xonepls) %*% Yonepls, type="heatmap")
for (idx_method in 1:length(methods_name)) {
  auto_t_path<-file.path("./simulation", methods_name[idx_method])
  if (!dir.exists(auto_t_path)){
    dir.create(auto_t_path)
  }
  if(methods_name[idx_method]=="PM_PLS"){
    train_model <- PLS_ZOO_AUTOT(Xonepm, 
                                 Yonepm, 
                                 E,method = methods_name[idx_method], 
                                 # lambda=c(0.05),
                                 lambda=lambda_set,
                                 k=5)
    # plot_lollipop(train_model,"PM_PLS","u")
    # plot_lollipop(train_model,"PM_PLS","v")
  }else if(methods_name[idx_method]=="PLS"){
    train_model <- PLS_ZOO_AUTOT(Xonepls, 
                                 Yonepls, 
                                 E,method = methods_name[idx_method], 
                                 k=5)
    # plot_lollipop(train_model,"PLS","u")
    # plot_lollipop(train_model,"PLS","v")
  }
  
  saveRDS(train_model,file = paste0(auto_t_path,"/best_train.rds"))
  saveRDS(X,file = paste0(auto_t_path,"/Simu_X_2.rds"))
  saveRDS(Y,file = paste0(auto_t_path,"/Simu_Y_2.rds"))
  saveRDS(E,file = paste0(auto_t_path,"/Simu_E_2.rds"))
  saveRDS(X_train_val,file = paste0(auto_t_path,"/Simu_X_train_val_2.rds"))
  saveRDS(Y_train_val,file = paste0(auto_t_path,"/Simu_Y_train_val_2.rds"))
  saveRDS(X_test,file = paste0(auto_t_path,"/Simu_X_test_2.rds"))
  saveRDS(Y_test,file = paste0(auto_t_path,"/Simu_Y_test_2.rds"))
  ## Testing performance
  res_table[1,idx_method]<-cor(X_test%*%train_model$bestModel$u,Y_test%*%train_model$bestModel$v)[1]
}

