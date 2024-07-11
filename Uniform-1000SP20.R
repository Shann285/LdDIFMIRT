rm(list=ls()) #clear screen

##set working dir
setwd("C:/Users/dell/Desktop")

library(MASS)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


bt_to<-date()
set.seed(12)
# specify the item design matrix and factor for model input
N <- 1000
P <- 15
D <- 2
nc <- 4
CNUM <- 50
YFILE <- "Y.txt"
XFILE <- "X.txt"

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
1,1),nrow = P,ncol = D,byr=T)
M <- sum(Ldes)

v0<-matrix(rnorm(P,0,1), nrow=P,ncol=1,byr=T)

K<-matrix(c(
0,0,0,0,
0,0,0,0,
0,0,0,0,
0.3,0,0,0,
0,0,0,0,
0,0,0,0,
0,0,0,0,
0,0.3,0,0,
0,0,0,0,
0,0,0,0,
0,0,0,0,
0,0,0,0,
0.3,0.3,0,0,
0,0,0,0,
0,0,0,0),nrow=P,ncol=nc,byr=T)

L0<-matrix(c(
1,0,
0,1,
1.4,0,
1.2,0,
0.8,0,
0.6,0,
0,1.4,
0,1.2,
0,0.8,
0,0.6,
0.6,1.4,
0.8,1.2,
1.0,1.0,
1.2,0.8,
1.4,0.6),nrow=P,ncol=D,byr=T)


for(CIR in 1:CNUM){
bb<-matrix(0, nrow=N, ncol=nc)
bb[,1]<- rnorm(N, 0, 1)                 ##rnorm(N, 0, 1)
bb[,2]<- rbinom(N, 1, 0.5)
bb[,3:4]<- mvrnorm(N, c(0,0), Sigma=matrix(c(1,0.5,0.5,1),nrow=2))

aa<-matrix(0, nrow=N, ncol=D)
for(n in 1:N){
   m1<- 0 + 0.5*bb[n,2]
   m2<- 0 - 0.5*bb[n,2]
   sig1<- 1*exp(0.0*bb[n,2])
   sig2<- 1*exp(0.0*bb[n,2])
   sigg<-diag(c(sig1,sig2))%*%matrix(c(1,0.5,0.5,1),nrow=D)%*%diag(c(sig1,sig2))
   aa[n,]<-mvrnorm(1, c(m1,m2), Sigma=sigg)
}


##aa<-mvrnorm(N, rep(0,D), Sigma=diag(1,D))

vi= matrix(rep(v0,N),ncol=N) + K%*%t(bb)

li= L0%*%t(aa) 

y=matrix(0,nrow=P,ncol=N)               
for(i in 1:N){
  for(j in 1:P){
    y[j,i]=rbinom(1, 1, pnorm(vi[j,i]+li[j,i]) )
  }
}

y=t(y)
x=bb

write(y, file="Y.txt", ncol=dim(y)[2], append=T)
write(x, file="X.txt", ncol=dim(x)[2], append=T)
print(CIR)
}


Anch = c(0,0,rep(1,P-2))

for(CIR in 1:CNUM){
 bt<-proc.time()
  y<-matrix(0,nrow=N,ncol=P)
  y<-matrix(scan(YFILE, skip=(CIR-1)*N, nlines=N), nrow=N, ncol=P)
  x<-matrix(0,nrow=N,ncol=nc)
  x<-matrix(scan(XFILE, skip=(CIR-1)*N, nlines=N), nrow=N, ncol=nc)

# model input and starting values
fa.data <- list(P =P, N =N, Y =y, D=D, nc =nc, Ldes = Ldes, M =M, Anch =Anch,
                    LNprior = 2,
                    X = x,
                    gamma_a = 9, gamma_b = 3)

init_model = function(){
  init.values <- list(
    Lr = rep(0.1, (M-2)) + runif((M-2),0,.1),
    nu = rep(0.1, P) + runif(P,0,.1),
    nu_dif = matrix((rep(0.1, (P-2)*nc) + runif((P-2)*nc,0,.1)), nrow=(P-2), ncol=nc),
    mu_imp = matrix((rep(0.1, D*nc) + runif(D*nc,0,.1)), nrow=D, ncol=nc),
    phi_imp = matrix((rep(0.1, D*nc) + runif(D*nc,0,.1)), nrow=D, ncol=nc),
    tau = rep(1, D),
    L_Omega = matrix(c(1,0,0,1),nrow=D))
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
cholesky_factor_corr[D] L_Omega;
}

transformed parameters {
vector[D] fac_scor[N]; // person factor scores
matrix[P,D] L; //  constrained factor loading matrix
matrix[P,nc] nu_1;
matrix[N,P] mu;
matrix[D,D] Sig_fac;
Sig_fac <- L_Omega * L_Omega';
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
  fac_scor[k] <- rep_vector(0,D) + mu_imp[,1]* X[k,1] + mu_imp[,2]* X[k,2] + mu_imp[,3]* X[k,3] + mu_imp[,4]* X[k,4] +
   diag_matrix(tau .*exp(phi_imp[,1]* X[k,1]+phi_imp[,2]* X[k,2]+phi_imp[,3]* X[k,3]+phi_imp[,4]* X[k,4])) * L_Omega * fac_dist_helper[k]';
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
L_Omega ~ lkj_corr_cholesky(2);
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
                                 "mu_imp", "phi_imp", "tau","Sig_fac"),
                         chains = 3, iter = 4000, init = init_model, cores = 3)


aa<-summary(stan_ssp, probs = c(0.025, 0.975), pars = c("Lr",
                                 "nu","nu_dif", 
                                 "lambda12", "lambda1_n", 
                                 "mu_imp", "phi_imp", "tau","Sig_fac"))

#a1<-aa$summary[,c(4)]*aa$summary[,c(5)]
anum<-get_num_divergent(stan_ssp)
write(anum, file="ndivergent.txt", ncol=1, append=TRUE, sep="\t")
write(aa$summary[,c(7)], file="Rhh.txt", ncol=length(aa$summary[,c(7)]), append=TRUE, sep="\t")
#write(a1, file="Sign.txt", ncol=length(a1), append=TRUE, sep="\t")
write(aa$summary[,c(1)], file="mean.txt", ncol=length(aa$summary[,c(1)]), append=TRUE, sep="\t")
write(aa$summary[,c(2)], file="se_mean.txt", ncol=length(aa$summary[,c(2)]), append=TRUE, sep="\t")
write(aa$summary[,c(3)], file="sd.txt", ncol=length(aa$summary[,c(3)]), append=TRUE, sep="\t")
write(aa$summary[,c(4)], file="quant1.txt", ncol=length(aa$summary[,c(4)]), append=TRUE, sep="\t")
write(aa$summary[,c(5)], file="quant2.txt", ncol=length(aa$summary[,c(5)]), append=TRUE, sep="\t")

print(CIR)
et<-proc.time()
print((et-bt)[3])
write((et-bt)[3], file="time.txt", ncol=1, append=TRUE, sep="\t")
}

date()
save.image(paste("BaLassoMIRT",".RData",sep=""))
