if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR")
BiocManager::install("Heatplus")
library("Heatplus")
library(CCA)
library(mvtnorm)
library(lattice)
library(rstan)
library(edgeR)
library(PMA)
library(invgamma)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%% Model contains the Stan Model
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


## We add S1 and S2 as sequencing or library size normalization factor in the model as 
## they can adjust the potential disproportional number of reads in different samples. 
## They can be calculated from calcNormFactors function in edge R package (Robiinson et al., 2010)

Modelhorse <- "
data {
int<lower=1> N; // sample size
int<lower=1> n1;// n1 features in dat1
int<lower=1> n2;// n1 features in dat1
int<lower=1> d; // lower dimens. 
row_vector<lower=0>[N] S1;
row_vector<lower=0>[N] S2;
int X[N,n1];    // dat1
int Y[N,n2];    // dat2
real<lower=0> scale_global ; // scale for the half -t prior for tau // ( tau0 = scale_global * sigma )
real<lower=1> nu_global;    // degrees of freedom for the half-t prior for tau
real<lower=1> nu_local;    // degrees of freedom for the half-t priors for lambdas // ( nu_local = 1 corresponds to the horseshoe )
}

parameters {
real<lower=0> sig_te; // sd errors
real<lower=0> sig_lb; // sd errors
matrix[N,d] Z;
matrix[N,n1] teta;
matrix[N,n2] lambd;
cholesky_factor_cov[n1,d] A; 
cholesky_factor_cov[n2,d] B; 

row_vector[n1] mux;
row_vector[n2] muy;

real<lower=0> taua[n1];              // global shrinkage parameter
matrix<lower=0>[n1,d] lamba ; // local shrinkage parameters

real<lower=0> taub[n2];             // global shrinkage parameter
matrix<lower=0>[n2,d] lambb; // local shrinkage parameters
}

// transformed param. parts

transformed parameters {

real<lower=0> sigte_sr;
real<lower=0> siglb_sr;
matrix[N,n1] mu_te;
matrix[N,n2] mu_lb;

mu_te = Z*A';
mu_lb = Z*B';

siglb_sr = sqrt(sig_lb);
sigte_sr = sqrt(sig_te);

}

// Model Part
model{

// half-t prior for tau
taua ~ cauchy(0, scale_global);
taub ~ cauchy(0, scale_global);

sig_lb ~ scaled_inv_chi_square(10, .05);
sig_te ~ scaled_inv_chi_square(10, .05);


mux ~ normal(0, 10);
muy ~ normal(0, 10);


for(i in 1:n1)
{lamba[i] ~ cauchy(0, 1);
A[i] ~ normal(0, lamba[i] * taua[i]);
}

for(i in 1:n2)
{lambb[i] ~ cauchy(0, 1);
B[i] ~ normal(0, lambb[i] * taub[i]);
}

for(i in 1:N)
{
  Z[i] ~ normal(0, 1);
  teta[i] ~ normal(mu_te[i] + mux, sigte_sr);
  X[i] ~ poisson_log(teta[i] + log(S1[i]));
  
  lambd[i] ~ normal(mu_lb[i] + muy, siglb_sr);
  Y[i] ~ poisson_log(lambd[i] + log(S2[i]));
}
}
"

CompModel <- stan_model(model_code=Modelhorse)

#******************************************************#
#******************************************************#
X =   ## Store the X data sets
Y =   ## Store the Y data sets
 
X_fac <- calcNormFactors(X) # S1
Y_fac <- calcNormFactors(Y) # S2


##########################################################################

Poste_summaries<-function(out,d1,n1,n2) ## take in the MCMC output
{
  
  A=matrix(out[1:(d1*n1)],nrow=n1,ncol=d1,byrow=FALSE)
  
  ul=n1*d1+(1:(n2*d1)) 
  
  B=matrix(out[ul],nrow=n2,ncol=d1,byrow=FALSE)
  
  siglb = out[d1*n1+n2*d1+2]
  
  sigte = out[d1*n1+n2*d1+1]
  
  out_all_log =as.vector(Post_cor_logn(A,B,sigte,siglb)) 
  
  return(out_all_log)    
  
} 

## calculate canonical correlation from the posterior 

Post_cor <- function(A,B,sigte,siglb){
  
  n1=nrow(A)
  
  n2=nrow(B) 
  
  R12 <- tcrossprod(B,A)
  R11 <- tcrossprod(A) + diag(sigte,n1)
  R22 <- tcrossprod(B) + diag(siglb,n2)
  
  R21 <- t(R12)
  
  e1 <- solve(R22) %*% R12 %*% solve(R11) %*% R21
  canon.corr <- sqrt(eigen(e1)$values)
  return(canon.corr)
}

## calculate the correlation from the posterior

Post_cor_logn <-function(A,B,sigte,siglb) 
{
  
  n1=nrow(A)
  
  n2=nrow(B)
  
  zxl=rbind(tcrossprod(A) + diag(sigte,n1), tcrossprod(B, A))
  
  zxu=rbind(tcrossprod(A,B) , tcrossprod(B) + diag(siglb,n2)) 
  
  mat_r <- cbind(zxl,zxu)
  
  ouf=cov2cor(mat_r)
  
  return(ouf)   
  
}

## calculate correlation of the raw data 

cor_XY <- function(X,Y) 
{
  zxl=cbind(cov(t(X)), cov(t(X),t(Y)))
  
  zxu=cbind(cov(t(Y),t(X)), cov(t(Y))) 
  
  mat_r <- rbind(zxl,zxu)
  
  ouf=cov2cor(mat_r)
  
  return(ouf)
}

## confidence interval 

confidence_interval <- function(vector, interval) {
  # Standard deviation of sample
  vec_sd <- sd(vector)
  # Sample size
  n <- length(vector)
  # Mean of sample
  vec_mean <- mean(vector)
  # Error according to t distribution
  error <- qt((interval + 1)/2, df = n - 1) * vec_sd / sqrt(n)
  # Confidence interval as a vector
  result <- c("lower" = vec_mean - error, "upper" = vec_mean + error)
  return(result)
}


## Generate count data

XY_Simulate <- function(N,d,n1,n2,A,B,mux,muy,sigte,siglb){
  
  Z = matrix(rnorm(N*d),nrow = d) 
  
  tet = mux + A%*%Z + matrix(rnorm(n=n1*N,sd=sigte),ncol=N)
  lamb = muy + B%*%Z + matrix(rnorm(n=n2*N,sd=siglb),ncol=N)
  
  X = apply(exp(tet), c(1,2),rpois,n=1) # rpois(n=1,lambda = 10*exp(tet))
  Y = apply(exp(lamb), c(1,2),rpois,n=1)
  
  return(list("X"=X,"Y"=Y))
}


Data <- list("X"=t(X),"Y"=t(Y),"n1"=n1,"n2"=n2,"N"=N,"d"=d1,"S1"=X_fac,"S2"=Y_fac,"scale_global"=1,"nu_global"=1,"nu_local"=1) 
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


fithorse = sampling(CompModel, data=Data, pars=c("A","B","sig_te","sig_lb","Z"),chains=2,iter=10000,control = list(adapt_delta = 0.99,max_treedepth = 15))

outval=as.matrix(fithorse)

out_cor = t(apply(X=outval,MARGIN=1,FUN=Poste_summaries,n1=n1,n2=n2,d=d1))

mean_frob=mean(unlist(frob))
CI_frob=confidence_interval(unlist(frob),0.95)

mean_stein=mean(unlist(stein))
CI_stein=confidence_interval(unlist(stein),0.95)

#X is the mean posterior correlation
#X_tr: is th true correlation:

#frob = sqrt(sum(((as.vector(x) - as.vector(x_tr))^2)))   # Frobenius loss
#vx <- svd(cor1[[1]])
#mt <- Tl_corpcan[[1]]%*%(vx$u%*%diag(1/vx$d)%*%t(vx$v))
#sten <- sum(diag(mt)) - log(det(mt)) - nrow(Tl_corpcan[[1]])   # Stein loss




