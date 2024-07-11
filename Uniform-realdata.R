rm(list=ls()) #clear screen

##set working dir
setwd("C:/Users/dell/Desktop")

library(MASS)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(foreign)


mydata=read.spss('Canada52.sav')

a<-matrix(unlist(mydata), ncol=55)
aa<-a[,-1]
aa[,-1]<-aa[,-1]-1
aa1<-aa[(aa[,1]==18)|(aa[,1]==19)|(aa[,1]==20)|(aa[,1]==21),] ##age=18 19 20 21
aa1<-na.omit(aa1)
aa1[,1]<- scale(aa1[,1])

b<-cbind(aa1[,1:2], aa1[,13], aa1[,32], aa1[,14:31], aa1[,33:54])



bt_to<-date()
set.seed(12)
# specify the item design matrix and factor for model input
N<-as.integer(843)		
P<-as.integer(42)		    
D <- 2
nc <- 2

Ldes <- matrix(c(
1,0,
0,1,
1,1,
1,1,
1,1,
1,1,
1,1,
1,1,
1,1,
1,1,
1,1,
1,1,
1,1,
1,1,
1,1,
1,1,
1,1,
1,1,
1,1,
1,1,
1,1,
1,1,
1,1,
1,1,
1,1,
1,1,
1,1,
1,1,
1,1,
1,1,
1,1,
1,1,
1,1,
1,1,
1,1,
1,1,
1,1,
1,1,
1,1,
1,1,
1,1,
1,1),nrow = P,ncol = D,byr=T)
M <- sum(Ldes)

Anch = c(0,0,rep(1,P-2))


bt<-proc.time()
  y<- b[, -c(1,2)]

  x<- b[, c(1,2)]


# model input and starting values
fa.data <- list(P =P, N =N, Y =y, D=D, nc =nc, Ldes = Ldes, M =M, Anch =Anch,
                    LNprior = 2,
                    X = x,
                    gamma_a = 27, gamma_b = 3)

init_model = function(){
  init.values <- list(
    Lr = rep(0.1, (M-2)) + runif((M-2),0,.1),
    nu = rep(0.1, P) + runif(P,0,.1),
    nu_dif = matrix((rep(0.1, (P-2)*nc) + runif((P-2)*nc,0,.1)), nrow=(P-2), ncol=nc),
    mu_imp = matrix((rep(0.1, D*nc) + runif(D*nc,0,.1)), nrow=D, ncol=nc),
    phi_imp = matrix((rep(0.1, D*nc) + runif(D*nc,0,.1)), nrow=D, ncol=nc),
    tau = rep(1, D),
    z01 = rep(0.1, (1+nc)) + runif((1+nc),0,.1) )
  return(init.values);
}

# compile the model
stan_m <- stan_model(model_code = "
data {
int<lower=1> N; // sample size
int<lower=1> P; // number of items
int Y[N,P]; // data matrix of order [N,P]
int<lower=1> D; // number of latent dimensions
int<lower=1> nc; // number of latent dimensions
matrix[P,D] Ldes; // design matrix for factor loadings
int<lower=1> M; // number of non-zero factor loadings
int Anch[P];
real<lower=0> LNprior; // non-DIF prior SD
matrix[N,nc] X;
real gamma_a;
real gamma_b;
}

parameters {
vector<lower=-0.5>[M-2] Lr; // non-zero factor loading
vector[P] nu; // intercept
matrix[(P-2),nc] nu_dif;
matrix[N,D] fac_dist_helper; // helper for non-centered sampling
matrix[D,nc] mu_imp; // factor mean impact
matrix[D,nc] phi_imp; // factor sd impact
vector<lower=0>[M-2] lambda12; // Laplace variance on loading DIF
matrix<lower=0>[(P-2),nc] lambda1_n; // Laplace variance on intercept DIF
vector<lower=0>[D] tau;
vector[(1+nc)] z01;  // covariance effects
}

transformed parameters {
vector[D] fac_scor[N]; // person factor scores
matrix[P,D] L; //  constrained factor loading matrix
matrix[P,nc] nu_1;
matrix[N,P] mu;
corr_matrix[D] L_Omega[N];
for(n in 1:N){
  L_Omega[n,1,1]<- 1;
  L_Omega[n,2,2]<- 1;
  L_Omega[n,1,2]<- (exp(2*(z01[1]+X[n]*z01[2:(1+nc)]))-1)/(exp(2*(z01[1]+X[n]*z01[2:(1+nc)]))+1);
  L_Omega[n,2,1]<- (exp(2*(z01[1]+X[n]*z01[2:(1+nc)]))-1)/(exp(2*(z01[1]+X[n]*z01[2:(1+nc)]))+1);}
L[1,1] <- 1;
L[1,2] <- 0;
L[2,1] <- 0;
L[2,2] <- 1;
L[3:P] <- rep_matrix(0,(P-2),D);
nu_1 <- rep_matrix(0,P,nc);


// assign factor loading and DIF param
{
int temp;
temp <- 0;
for(i in 1:D){
  for(j in 3:P){
    if(Ldes[j,i]==1){
      temp <- temp + 1;
      L[j,i] <- Lr[temp];}
  }
}
}
{
int temp;
temp <- 0;
  for(j in 1:P){
     if(Anch[j]==1){
        temp <- temp + 1;
        nu_1[j] <-  nu_dif[temp];}
  }
}
      
 
// build non - centered factor score multi - norm dist . for performance ;
// Stan needs to work with the cholesky factor of phi.
for(k in 1:N){
  fac_scor[k] <- rep_vector(0,D) + mu_imp[,1]* X[k,1] + mu_imp[,2]* X[k,2] + 
   diag_matrix(tau .*exp(phi_imp[,1]* X[k,1]+phi_imp[,2]* X[k,2])) * cholesky_decompose(L_Omega[k]) * fac_dist_helper[k]';
  for(l in 1:P){
    mu[k,l] <- nu[l] + nu_1[l]* X[k]' + L[l] * fac_scor[k];}
}
}

model {
// the priors
lambda12 ~ gamma(gamma_a, gamma_b);
to_vector(lambda1_n) ~ gamma(gamma_a, gamma_b);
Lr ~ double_exponential(0, 1 ./ sqrt(lambda12));
to_vector(fac_dist_helper) ~ normal(0, 1);
nu ~ normal(0, LNprior);
to_vector(nu_dif) ~ double_exponential(0, 1 ./ sqrt(to_vector(lambda1_n)));
to_vector(mu_imp) ~ normal(0, LNprior);
to_vector(phi_imp) ~ normal(0, LNprior);
tau ~ cauchy(0, 2.5);
z01 ~ normal(0, LNprior);
// The likelihood
for(j in 1:P){
  Y[ ,j] ~ bernoulli(Phi(mu[ ,j]));
}
}
", verbose = T )


stan_ssp <- sampling(stan_m,
                         data = fa.data,
                         pars =c("Lr",
                                 "nu","nu_dif", 
                                 "lambda12", "lambda1_n", 
                                 "mu_imp", "phi_imp", "tau","z01"),
                         chains = 3, iter = 5000, init = init_model, cores = 3)


et<-proc.time()
print((et-bt)[3])

aa<-summary(stan_ssp, probs = c(0.025, 0.975), pars = c("Lr",
                                 "nu","nu_dif", 
                                 "lambda12", "lambda1_n", 
                                 "mu_imp", "phi_imp", "tau","z01"))

