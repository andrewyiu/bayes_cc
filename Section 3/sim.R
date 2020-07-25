###################################
## Simulation study in section 3 ##
###################################
library(foreach)
library(doParallel)
library(CVST)
library(MASS)
library(doRNG)
library(survival)

n <- 2000
beta0 <- 0.3
eta <- 0.01
nu <- 2.0
cores <- 20
trials <- 2000
cl <- makeCluster(cores)
registerDoParallel(cl)

start_time <- Sys.time()
estimators <-foreach(k=1:trials, .packages="survival",
               .combine='rbind',.options.RNG=1) %dorng%{
  
  x <- rnorm(n) # generate x
  
  # weibull hazard to generate failure times
  u <- runif(n)
  T <- (-log(u)/(eta*exp(beta0*x)))^(1/nu)
  
  # censoring times
  surv <- rbinom(n, 1, 0.2) # 0.2 prob that individual survives to end (t=3)
  C <- numeric(n)
  for (i in 1:n) {
    if (surv[i] == 1) {
      C[i] <- 3
    } else {
      C[i] <- runif(1, 0, 3)
    }
  }
  
  delta <- as.numeric(T < C) # case indicators
  Y <- pmin(T,C) # observed survival times
  
  samp <- sample(1:n, size = n/25) # sampling the subcohort
  fullcases <- which(delta == 1) # indices of cases
  numsubcases <- length(which(delta[samp]==1)) # number of subcohort cases
  subco <- union(samp, fullcases) # "extended" subcohort: union of subcohort and remaining cases
  size_sub <- length(samp)
  
  # at-risk indicators
  risk <- matrix(0, nrow = length(fullcases), ncol = n)
  for (i in 1:length(fullcases)) {
    for (j in 1:n) {
      if (Y[j] >= Y[fullcases[i]]) {
        risk[i,j] <- 1
      }
    }
  }
  
  Ysub <- Y[subco]
  deltasub <- delta[subco]
  xsub <- x[subco]
  
  testfull <- list(time= as.numeric(Y),
                   status=delta,
                   x=x)
  #chen-lo weights
  cl_wt <- numeric(length(subco))
  for (k in 1:length(subco)) {
    if (deltasub[k] == 1) {
      cl_wt[k] <- 1
    } else {
      #number of controls divided by number of controls in subcohort
      cl_wt[k] <- (n-sum(delta))/(size_sub-numsubcases)
    }
  }
  testcl <- list(time= as.numeric(Ysub),
                 status=deltasub,
                 x=xsub)
  #kalbfleisch-lawless weights
  kl_wt <- numeric(length(subco))
  for (k in 1:length(subco)) {
    if (deltasub[k] == 1) {
      kl_wt[k] <- 1
    } else {
      kl_wt[k] <- 25
    }
  }
  testkl <- list(time= as.numeric(Ysub),
                 status=deltasub,
                 x=xsub)
  
  #prentice
  testpren <- list(time= as.numeric(Ysub),
                   status=deltasub,
                   x=xsub)
  sam <- as.numeric(subco %in% samp)
  
  # Metropolis-Hastings
  propvar <- (2^2)*coxph(Surv(time,status) ~ x,data=testcl, weights <- cl_wt)$var #proposal variance set to 4 times chen-lo estimate variance
  warmup <- 1000
  iter <- 20000
  
  current <- beta0 #initial value set to beta0
  
  storesub <- exp(xsub*current) #pre-calculate exponential factors for observed covariates
  partprop <- matrix(, nrow = n, ncol = 1)
  
  dir_wts <- rexp(length(subco), 1) 
  dir_wts <- dir_wts/sum(dir_wts) # sample dirichlet weights for bayesian bootstrap
  
  # fill in missing covariates
  fill_cov <- numeric(n) # assign indices of covariates to each individual in the cohort
  fill_cov[subco] <- 1:length(subco) # individuals with observed covariates 
  # individuals with missing covariates are assigned covariates values by resampling according to dirichlet weights
  fill_cov[-subco] <- sample(x=1:length(subco), size = (n-length(subco)), prob = dir_wts, replace = TRUE)
  store <- storesub[fill_cov]
  partprop[,1] <- store
                 
  denom <- risk %*% partprop # compute denominators of the factors in the cox likelihood
  
  hcurr <- prod(partprop[fullcases,1]*n/denom[,1]) # compute likelihood
  
  # warm-up iterations
  warm <- numeric(warmup)
  for (w in 1:warmup) {
    prop <- rnorm(n=1, mean =current, sd = sqrt(propvar)) # generate proposal for beta
    # compute proposal likelihood term
    storesub <- exp(xsub*prop)
    dir_wts <- rexp(length(subco), 1)
    dir_wts <- dir_wts/sum(dir_wts)
    fill_cov <- numeric(n)
    fill_cov[subco] <- 1:length(subco)
    fill_cov[-subco] <- sample(x=1:length(subco), size = (n-length(subco)), prob = dir_wts, replace = TRUE)
    store <- storesub[fill_cov]
    partprop[,1] <- store
    denom <- risk %*% partprop
    hprop <- prod(partprop[fullcases,1]*n/denom[,1])
    
    if (hprop > 0) {
      ratio <- hprop/hcurr # acceptance ratio
      a <- rbinom(n=1, 1, min(ratio,1))
      if (a == 1) {
        current <- prop
        hcurr <- hprop
      }
    }
    warm[w] <- current
  }
  post_samp <- numeric(iter) # posterior samples
  for (t in 1:iter) {
    prop <- rnorm(n=1, mean =current, sd = sqrt(propvar))
    storesub <- exp(xsub*prop)
    dir_wts <- rexp(length(subco), 1)
    dir_wts <- dir_wts/sum(dir_wts)
    fill_cov <- numeric(n)
    fill_cov[subco] <- 1:length(subco)
    fill_cov[-subco] <- sample(x=1:length(subco), size = (n-length(subco)), prob = dir_wts, replace = TRUE)
    store <- storesub[fill_cov]
    partprop[,1] <- store
    denom <- risk %*% partprop
    
    hprop <- prod(partprop[fullcases,1]*n/denom[,1])
    
    if (hprop > 0) {
      ratio <- hprop/hcurr
      a <- rbinom(n=1, 1, min(ratio,1))
      if (a == 1) {
        current <- prop
        hcurr <- hprop
      }
    }
    post_samp[t] <- current
  }
  # check whether intervals contain the truth
  samp_ord <- post_samp[order(post_samp)]
  if (beta0 > samp_ord[(iter*0.025)] && beta0 <= samp_ord[(iter*0.975)]) {
    cover1 <- 1
  } else {
    cover1 <- 0
  }
  # coverage for full data
  fit0 <- coxph(Surv(time,status) ~ x,data=testfull)
  est0 <- fit0$coefficients
  if (beta0 > (est0-qnorm(0.975)*sqrt(fit0$var)) && beta0 <= (est0+qnorm(0.975)*sqrt(fit0$var))) {
    cover0 <- 1
  } else {
    cover0 <- 0
  }
  
  est1 <- mean(post_samp)
  # coverage for chen-lo
  fit2 <- coxph(Surv(time,status) ~ x,data=testcl, weights <- cl_wt, robust = TRUE)
  est2 <- fit2$coefficients
  if (beta0 > (est2-qnorm(0.975)*sqrt(fit2$var)) && beta0 <= (est2+qnorm(0.975)*sqrt(fit2$var))) {
    cover2 <- 1
  } else {
    cover2 <- 0
  }
  # coverage for kalbfleisch-lawless
  fit3 <- coxph(Surv(time,status) ~ x,data=testkl, weights <- kl_wt, robust = TRUE)
  est3 <- fit3$coefficients
  if (beta0 > (est3-qnorm(0.975)*sqrt(fit3$var)) && beta0 <= (est3+qnorm(0.975)*sqrt(fit3$var))) {
    cover3 <- 1
  } else {
    cover3 <- 0
  }
  #coverage for prentice
  fit4 <- cch(Surv(time,status) ~ x,data=testpren, subco = sam, id = 1:length(subco), cohort.size=length(subco), method = "Prentice")
  est4 <- fit4$coefficients
  if (beta0 > (est4-qnorm(0.975)*sqrt(fit4$var)) && beta0 <= (est4+qnorm(0.975)*sqrt(fit4$var))) {
    cover4 <- 1
  } else {
    cover4 <- 0
  }
  
  c(est0,est1,est2,est3,est4, cover0, cover1, cover2, cover3, cover4)
}
end <- Sys.time()-start_time
saveRDS(estimators, file = "estimators.rds")
print(end)
print(c(beta0, eta, nu))

# printing table 1 values
coverage <- colMeans(estimators[,6:10]*100)
values <- estimators[,1:5]
bias <- colMeans(values)-beta0
stan_dev <- apply(values,2,sd)
mse <- bias^2+stan_dev^2
rmse <- round(sqrt(mse),3)
bias <- round(bias, 3)
stan_dev <- round(stan_dev,3)
re <- mse[1]/mse
names <- c("Full", "Bayes", "CL", "KL", "Prentice")
tab <- data.frame(cbind(names), bias, stan_dev, rmse, re, coverage)
print(tab) 

stopCluster(cl)


