rm(list=ls()) #clear screen

##set working dir
setwd("C:/Users/dell/Desktop")

library(MASS)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
# ---- Necessary packages for EML1----
library(glmnet)     # L1-logistic regression
library(mvtnorm)    # dmvnorm
library(progress)   # progress bar


bt_to<-date()
set.seed(12)
# specify the item design matrix and factor for model input
N <- 500
P <- 10
D <- 2
nc <- 1
CNUM <- 50
YFILE <- "Y.txt"
XFILE <- "X.txt"

v0<-matrix(0, nrow=P,ncol=1,byr=T)

K<-matrix(c(
0,
0,
0,
0.3,
0,
0,
0,
0.3,
0,
0),nrow=P,ncol=nc,byr=T)

L0<-matrix(c(
1,0,
0,1,
1,0,
1,0,
1,0.3,
1,0.3,
0,1,
0,1,
0.3,1,
0.3,1),nrow=P,ncol=D,byr=T)



for(CIR in 1:CNUM){
bb<-matrix(0, nrow=N, ncol=nc)
bb[,1]<- c(rep(0,N/2), rep(1,N/2))              ##rnorm(N, 0, 1)


aa<-matrix(0, nrow=N, ncol=D)
for(n in 1:N){
   m1<- 0 + 0.5*bb[n,1]
   m2<- 0 - 0.5*bb[n,1]
   sig1<- 1*exp(0.0*bb[n,1])
   sig2<- 1*exp(0.0*bb[n,1])
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

##EML1

# ---- Some useful functions ----
Grid_pts <- function(K, lb=-4, ub=4, np=11){
  # K     : no. of latent variables.
  # lb, ub: lower bound, upper bound.
  # np    : no. of grid points for each variable.
  # output: grid points.
  
  p_list <- list()
  for(k in 1:K){
    p_list[[k]] <- seq(from=lb, to=ub, length.out=np)
  }
  grid_pts <- as.matrix(expand.grid(p_list, stringsAsFactors=FALSE))
  colnames(grid_pts) <- NULL
  return(grid_pts)
}


Compute_p_tilde <- function(y, x_grid, A, b, mu, sigma){
  # y        : responses of all subjects, a matrix with N*J.
  # x_grid   : grid points of x, note that all subjects have same grid points.
  # A, b     : current item parameters.
  # mu, sigma: parameters of p.d.f. of x grid.
  # return   : p tilde.
  
  M <- nrow(x_grid)
  N <- nrow(y)
  
  log_px  <- dmvnorm(x_grid, mean=mu, sigma=sigma, log=T)
  logit_p <- x_grid%*%A + rep(1, M)%*%t(b)
  p       <- plogis(logit_p)
  
  prob <- rep(0, N*M)
  for(i in 1:N){
    log_pyx <- rowSums(dbinom(rep(1, M)%*%t(y[i, ]), 1, p, log=T))
    loglik  <- log_px + log_pyx
    lik     <- exp(loglik)
    prob[1:M+(i-1)*M] <- lik/sum(lik)
  } 
  return(prob)
}


# -------- EML1 M2PL function --------
EM_MIRT <- function(y,            # data set, all responses of all subjects
                    mu_t,         # true mu (mean vector of latent traits)
                    sigma_t,      # true sigma (covariance of latent traits)
                    fixed,        # which item not do sub-model selection
                    lambda_list,  # penalized parameters list
                    A_init,       # initial value of A for EM algorithm
                    b_init,       # initial value of b for EM algorithm
                    grid_num      # number of grid points for each latent trait
){
  
  # ---- Dimension of A & y ----
  J <- ncol(A_init)  # number of items
  K <- nrow(A_init)  # number of latent traits
  N <- nrow(y)       # number of subjects
  
  # ---- Grid points ----
  x_grid <- Grid_pts(K=K, lb=-4, ub=4, np=grid_num)
  M      <- nrow(x_grid)
  
  y_samp <- apply(y, MARGIN=2, FUN=function(x){rep(x,each=M,times=1)})
  p_samp <- rep(NA, times=N*M)
  x_samp <- matrix(0,N*M,K)    
  for(i in 1:N){
    x_samp[(1:M)+(i-1)*M,] <- x_grid
  }
  
  # ---- Simulation settings ----
  nla    <- length(lambda_list)
  b_list <- list()
  A_list <- list()
  for(ila in 1:nla){
    b_list[[ila]] <- rep(0, J)
    A_list[[ila]] <- matrix(data=0, nrow=K, ncol=J)
  }
  bic_vec      <- rep(Inf, nla)
  cr_vec       <- rep(0, nla)
  time_em1_vec <- rep(0, nla)
  time_em2_vec <- rep(0, nla)
  iter_em1_vec <- rep(0, nla)
  iter_em2_vec <- rep(0, nla)
  
  bics <- rep(Inf, J)
  
  # ---- Simulation ----
  # count <- 0
  time_total <- proc.time()
  
  for(ila in 1:nla){
    
    cat(sprintf("---- ila: %02d ----\n", ila))
    
    # ---- EM L1 initialization ----
    if(ila==1){
      A_c   <- A_init
      Mod_c <- A_c!=0
      b_c   <- b_init
      
      A_new   <- A_c
      Mod_new <- Mod_c
      b_new   <- b_c
    }
    
    # ---- EM L1 iterations ----
    t_em1 <- proc.time()
    iter_em1 <- 0
    
    while(iter_em1 < 100){
      
      iter_em1 <- iter_em1 + 1
      pb_em <- progress_bar$new(format=sprintf("la:%02d EM-L1:%03d [:bar] :percent eta::eta", ila, iter_em1),
                                total=J+1, clear=TRUE, width=60, show_after=0)
      pb_em$tick(0)  # progress bar
      
      # ---- E-step: ----
      p_samp <- Compute_p_tilde(y, x_grid, A_c, b_c, mu_t, sigma_t)
      pb_em$tick(1)
      
      # ---- M-step: ----
      for(j in fixed){
        
        excl <- which(A_init[,j]==0)
        
        fit <- glmnet(x=x_samp, y=y_samp[,j], weights=p_samp, family="binomial",
                      alpha=1, lambda=0, standardize=FALSE, exclude=excl)
        
        bics_j    <- (1-fit$dev.ratio)*fit$nulldev + log(N)*fit$df
        A_new[,j] <- as.vector(coef(fit)[-1])
        b_new[j]  <- as.numeric(coef(fit)[1])
        bics[j]   <- bics_j
        
        pb_em$tick(1)
      }
      
      for(j in (1:J)[-fixed]){
        
        fit <- glmnet(x=x_samp, y=y_samp[,j], weights=p_samp, family="binomial",
                      alpha=1, lambda=lambda_list[ila], standardize=FALSE)
        
        bics_j    <- (1-fit$dev.ratio)*fit$nulldev + log(N)*fit$df
        A_new[,j] <- as.vector(coef(fit)[-1])
        b_new[j]  <- as.numeric(coef(fit)[1])
        bics[j]   <- bics_j
        pb_em$tick(1)
      }
      
      Mod_new <- A_new!=0
      # ---- Display the new parameters ----
      cat("A_new:\n");   print(A_new)
      
      # ---- Stop criterion ----
      if(all(Mod_new==Mod_c)){
        err <- max(abs(A_new[Mod_new]-A_c[Mod_c])/abs(A_c[Mod_c])) # maximum relative difference
        cat("err:", err, "\n")
        if(err < 0.1){
          break
        }
      }
      else{
        cat("Models wasn't same.\n")
      }
      
      # ---- Replace current parameter ----
      A_c   <- A_new
      b_c   <- b_new
      Mod_c <- Mod_new
      
    } # end while EM-L1
    
    t_em1 <- proc.time() - t_em1
    cat("t_em1:", as.numeric(t_em1[3]), "seconds.\n")
    
    time_em1_vec[ila] <- as.numeric(t_em1[3])
    iter_em1_vec[ila] <- iter_em1
    
    
    # ---- EM parameters estimate ----
    # ---- Note that: if the model is same to previous lambda, skip this one. 
    if(ila > 1){
      if(all(Mod_new==(A_list[[ila-1]]!=0))){
        
        A_list[[ila]] <- A_list[[ila-1]]
        b_list[[ila]] <- b_list[[ila-1]]
        bic_vec[ila]  <- bic_vec[ila-1]
        next
      }
    }
    
    A_c <- A_new
    b_c <- b_new
    
    t_em2 <- proc.time()
    iter_em2 <- 0
    while(iter_em2 < 100){
      
      iter_em2 <- iter_em2 + 1
      pb_em <- progress_bar$new(format=sprintf("la:%02d EM:%03d [:bar] :percent eta::eta", ila, iter_em2),
                                total=J+1, clear=TRUE, width=60, show_after=0)  
      pb_em$tick(0)  # progress bar
      
      # ---- E-step: ----
      p_samp <- Compute_p_tilde(y=y, x=x_grid, A=A_c, b=b_c, mu=mu_t, sigma=sigma_t)
      pb_em$tick(1)
      
      # ---- M-step: ----
      for(j in 1:J){
        
        if(sum(A_c[,j]==0)==K){
          # -- if jth column are all zeros, skip this estimation --
          A_new[,j] <- A_c[,j]
          b_new[j]  <- b_c[j]
        }
        else{
          excl <- which(A_c[,j]==0)
          fit <- glmnet(x=x_samp, y=y_samp[,j], weights=p_samp, family="binomial",
                        alpha=1, lambda=0,  standardize=FALSE, exclude=excl)
          
          bics_j    <- (1-fit$dev.ratio)*fit$nulldev + log(N)*fit$df
          A_new[,j] <- as.vector(coef(fit)[-1])
          b_new[j]  <- as.numeric(coef(fit)[1])
          bics[j]   <- bics_j
        }
        
        pb_em$tick(1)
      }
      
      # ---- Display the new parameters ----
      cat("A_new:\n");   print(A_new)
      
      # ---- Stop criterion ----
      # -- Note: Why we check whether the models same?
      #          One time, the non-zero parameter turn to zero through estimate.
      
      if(all((A_new!=0)==(A_c!=0))){
        err <- max(abs(A_new[A_new!=0]-A_c[A_c!=0])/abs(A_c[A_c!=0])) # maximum relative difference
        cat("err:", err, "\n")
        if(err < 0.1){
          break
        }
      }
      else{
        cat("Models wasn't same.\n")
      }
      
      # ---- Replace current parameter ----
      A_c <- A_new
      b_c <- b_new
      
    } # end while EM
    
    t_em2 <- proc.time() - t_em2
    cat("em2_use:", as.numeric(t_em2[3]), "seconds.\n")
    
    time_em2_vec[ila] <- as.numeric(t_em2[3])
    iter_em2_vec[ila] <- iter_em2
    
    
    # ---- Calculate bic ----
    # bic <- Compute_Bic(y, x_grid, A_new, b_new, prob_samp)
    bic <- sum(bics)
    
    # ---- Collect result ----
    A_list[[ila]] <- A_new
    b_list[[ila]] <- b_new
    bic_vec[ila]  <- bic
    
    A_t_sub       <- A_t[,-fixed]
    A_opt_sub     <- A_list[[ila]][,-fixed]
    cr_vec[ila]   <- sum((A_t_sub!=0)==(A_opt_sub!=0))/(J-K)/K
    
  } # end for(ila in 1:nla)
  
  
  # ---- Model selection ----
  time_total <- proc.time() - time_total
  time_total <- as.numeric(time_total[3])
  cat("use_total:", time_total, "\n")
  
  bic_min <- which.min(bic_vec)
  A_opt <- A_list[[bic_min]]
  b_opt <- b_list[[bic_min]]
  
  # ---- Return output ----
  result <- list(A_list=A_list,
                 b_list=b_list,
                 bic_vec=bic_vec,
                 bic_min=bic_min,
                 time_em1_vec=time_em1_vec,
                 time_em2_vec=time_em2_vec,
                 iter_em1_vec=iter_em1_vec,
                 iter_em2_vec=iter_em2_vec,
                 time_total=time_total,
                 A_opt=A_opt,
                 b_opt=b_opt
  )
  return(result)
  
} # end EML1 MIRT function




Anch = c(0,0,rep(1,P-2))

for(CIR in 1:CNUM){
 bt<-proc.time()
  y<-matrix(0,nrow=N,ncol=P)
  y<-matrix(scan(YFILE, skip=(CIR-1)*N, nlines=N), nrow=N, ncol=P)
  x<-matrix(0,nrow=N,ncol=nc)
  x<-matrix(scan(XFILE, skip=(CIR-1)*N, nlines=N), nrow=N, ncol=nc)

  # ---- EM initial parameters ----
  A_t <- c(1,0,1,1,1,1,0,0,0.3,0.3,
           0,1,0,0,0.3,0.3,1,1,1,1)*1.7
  A_t <- matrix(A_t, nrow=2, ncol=10, byrow=T)
  A_init <- matrix(data=1/P, nrow=D, ncol=P, byrow=TRUE)
  fixed <- c(1,2)
  
  A_init[,fixed] <- matrix(c(1.7,0,0,1.7),nrow=D)
  b_init <- rep(0,P)
  mu_t <- rep(0,D)
  sigma_t <- matrix(0.5,D,D); diag(sigma_t) <- 1
  grid_num <- 11
  lambda_list <- seq(from=0.01, to=0.05, by=0.001)
  yref <- y
  # ---- EM L1 ----
  a<-EM_MIRT(y=yref,               
             mu_t=mu_t, 
             sigma_t=sigma_t,
             fixed=fixed,
             lambda_list,
             A_init=A_init,
             b_init=b_init,
             grid_num=grid_num
  )
  
  write(as.vector(a$A_opt), file="Aopt.txt", ncol=(P*D), append=TRUE, sep="\t")
  
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
    1,1),nrow = P,ncol = D,byr=T)
  
  Ldes[t(a$A_opt)==0]=0
  M <- sum(Ldes)
  
  
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
  for(j in 1:P){
    if(Ldes[j,i]==1 && Anch[j]==1){
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
  fac_scor[k] <- rep_vector(0,D) + mu_imp[,1]* X[k,1] +
   diag_matrix(tau .*exp(phi_imp[,1]* X[k,1])) * L_Omega * fac_dist_helper[k]';
  for(l in 1:P){
    mu[k,l] <- nu[l] + nu_1[l]* X[k]' + L[l] * fac_scor[k];}
}
}

model {
// the priors
to_vector(lambda1_n) ~ gamma(gamma_a, gamma_b);
Lr ~ normal(0, LNprior);
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
                         pars =c("L",
                                 "nu","nu_dif", 
                                 "lambda1_n",
                                 "mu_imp", "phi_imp", "tau","Sig_fac"),
                         chains = 3, iter = 4000, init = init_model, cores = 3)

aa<-summary(stan_ssp, probs = c(0.025, 0.975), pars = c("L",
                                                        "nu","nu_dif", 
                                                        "lambda1_n",
                                                        "mu_imp", "phi_imp", "tau","Sig_fac"))


a1<-aa$summary[,c(4)]*aa$summary[,c(5)]
anum<-get_num_divergent(stan_ssp)
write(anum, file="ndivergent.txt", ncol=1, append=TRUE, sep="\t")
write(aa$summary[,c(7)], file="Rhh.txt", ncol=length(aa$summary[,c(7)]), append=TRUE, sep="\t")
write(a1, file="Sign.txt", ncol=length(a1), append=TRUE, sep="\t")
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

