library(pracma)
#install.packages("R.matlab")
library(R.matlab)

n = 10; 
p = 8; 
q = 15;
set.seed(10)
z = c(2.9298929,  0.1160284 ,-0.6535712, -0.9627428,  0.9606273 ,-0.2462755 ,-0.1951486 ,0.9684572  ,1.6764075 , 0.6969553)

set.seed(11)
X = as.matrix(data.frame(z1=z, z2=z, z3=z, matrix(rnorm(n*(p-3)),nrow = n,ncol = p-3)))
X = readMat("D:/Box/02-Research/210101-CCA/210704-R-Code/X.mat")[[1]]

set.seed(12)
Y = as.matrix(data.frame(z1=z, z2=z, z3=z, z4=z, z5=z, matrix(rnorm(n*(q-5)),nrow = n,ncol = q-5)))
Y = readMat("D:/Box/02-Research/210101-CCA/210704-R-Code/Y.mat")[[1]]

E = matrix(0,p,q)
E[1:3,1:5] = 1

temp = PM_SPLS(X, Y, E);
temp = PM_PLS(X, Y, E);
temp = SPLS(X, Y);
temp = PLS(X, Y);

temp = PM_SPLS(X, Y, E, standardize="Raw");

temp = PM_SPLS(X, Y, E, c1=sqrt(3),c2=sqrt(5),lambda=0.3)
5
temp = PM_SPLS(X, Y, E, c1=sqrt(3),c2=sqrt(5),lambda=0.3,u0=matrix(rnorm(p),ncol = 1),v0=matrix(rnorm(q),ncol = 1))

## Test the autotune algorithm
temp <- PM_SPLS_AUTOT(X, Y, E, c1=c(1,2,3,4),c2=c(1,2,3,4),lambda=c(0.25,0.5,0.75),k=5)
temp <- PM_SPLS_AUTOT(X, Y, E, c1=c(1,2,3,4),c2=c(1,2,3,4),lambda=c(0.25,0.5,0.75),k=5,stratified = c(0,0,0,0,0,1,1,1,1,1))
temp <- PM_SPLS_AUTOT(X, Y, E,lambda=c(0.25,0.5,0.75),k=5)
temp <- PM_SPLS_AUTOT(X, Y, E,lambda=c(0.25,0.75),k=5, standardize="Raw")


standardize="Center"
c1=NULL
c2=NULL
lambda=0.25
u0=NULL
v0=NULL
maxiter=100
tolvar=10e-3
tolfun=10e-6
QUIET=F