###############################################################
## Time the section 3 simulations for different sample sizes ##
###############################################################
library(foreach)
library(doParallel)
library(CVST)
library(MASS)
library(doRNG)
library(survival)

beta0 <- 0.3
eta <- 0.01
nu <- 2.0
cores <- 20
trials <- 100
cl <- makeCluster(cores)
registerDoParallel(cl)
results <- matrix(,nrow = 10, ncol = trials)
for (m in 1:10) {
  n <- m*1000
times <-foreach(k=1:trials, .packages="survival",
                     .combine='rbind',.options.RNG=1) %dorng%{
                       start_time <- Sys.time()
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

sam <- as.numeric(subco %in% samp)
B <- 1

# chen-lo weights to use for proposal variance
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

end <- Sys.time()-start_time
c(end)


                     }
results[m,] <- times
}

saveRDS(results, file = "times.rds")

stopCluster(cl)

## plot figure 1
# timevec <- as.vector(t(times))
# mult <- numeric(0)
# for (m in 1:10) {
#   mult <- c(mult, rep(m,100))
# }
# timevalues <- seq(1000,10000,200)
# mult2 <- mult^2
# quad <- lm(timevec ~ mult + mult2)
# predictedcounts <- predict(quad,list(mult=timevalues, mult2=timevalues^2))
# plot(mult, timevec, xlab = "n", ylab = "Duration (s)")
# lines(timevalues, predictedcounts)
