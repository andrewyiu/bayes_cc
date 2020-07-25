###########################################
## Application to the EPIC-Norfolk study ##
###########################################

# load data
data = read.csv("mydata_up.csv")

for (j in 1:25653) {
  if (is.na(data$dmstatus_ver_outc[j]) == TRUE) {
    data$dmstatus_ver_outc[j] <- 3
  }
}

# load initial values for beta and estimate of posterior variance
prop_list <- readRDS("prop_list.rds")
init_est <- prop_list$init_est
V_est <- prop_list$V_est
library(CVST)
library(MASS)
library(survival)
library(MCMCpack)
library(MBSP)
set.seed(123)

## seeds for the trace plots
#set.seed(234)
#set.seed(345)

# correlation parameters
rho <- 0.995
rhoxi <- 0.995

# total number of MCMC iterations
iter <- 100

start_time <- Sys.time()

# remove prevalent cases
mydat <- data[-which(data$dmstatus_ver == 2),]

# remove individuals with unknown status
mydat <- mydat[-which(mydat$dmstatus_ver == -1),]

# remove individuals with missing potential confounders
mydat <- mydat[which(is.na(mydat$pa_index)==FALSE),]
mydat <- mydat[which(is.na(mydat$bmi_adj)==FALSE),]
mydat <- mydat[which(is.na(mydat$waist_adj)==FALSE),]
mydat <- mydat[which(is.na(mydat$sex)==FALSE),]
mydat <- mydat[which(is.na(mydat$age_recr)==FALSE),]

# remove excleier variables
mydat <- mydat[-which(is.na(mydat$excleier)==TRUE),]
mydat <- mydat[-which(mydat$excleier == 1),]
mydat <- mydat[-which(mydat$excleier == 3),]

n <- dim(mydat)[1]
Y <- mydat$fup_time
delta <- as.integer(mydat$dmstatus_ver_outc == 1 | mydat$dmstatus_ver_outc == 2)

# indices of the individuals who are incident cases
fullcases <- which(mydat$dmstatus_ver_outc == 1 | mydat$dmstatus_ver_outc == 2)

# individuals with fatty acid measurements
subco <- which(is.na(mydat$C100)==FALSE)

num_miss <- n-length(subco)
sam <- which(mydat$dmstatus_ver_outc == 1 | mydat$dmstatus_ver_outc == 0) #individuals sampled into subcohort

z <- matrix(, nrow = dim(mydat)[1], ncol = 9) # saturated fatty acids
unsat <- 100-rowSums(mydat[,c(35, 39, 33, 37, 41, 47, 53, 56, 59)]) # total of unsaturated fatty acids
z[, 1] <- mydat[,c(35)] # C15:0
z[, 2] <- mydat[,c(39)] # C17:0
z[, 3] <- mydat[,c(33)] # C14:0
z[, 4] <- mydat[,c(37)] # C16:0
z[, 5] <- mydat[,c(41)] # C18:0
z[, 6] <- mydat[,c(47)] # C20:0
z[, 7] <- mydat[,c(53)] # C22:0
z[, 8] <- mydat[,c(56)] # C23:0
z[, 9] <- mydat[,c(59)] # C24:0

w <- matrix(, nrow = dim(mydat)[1], ncol = 7) # potential confounders
w[,1] <- mydat[,c(13)]-1 # sex
w[,2] <- mydat[,c(15)] # age at recruitment
w[,3] <- mydat[,c(16)] # adjusted waist
w[,4] <- mydat[,c(17)] # adjusted bmi

# physical activity index
w[,5] <- as.numeric(mydat[,c(19)]==1) # "inactive"
w[,6] <- as.numeric(mydat[,c(19)]==2) # "moderately inactive"
w[,7] <- as.numeric(mydat[,c(19)]==3) # "moderately active"

sdw <- apply(w,2,sd)
w <- w %*% diag((1/c(1,sdw[2:4],1,1,1))) # normalize by sd's


# additive logratio transformation (replace any zero fatty acid measurements with 0.005 first)
for (i in 1:length(subco)) {
  for (j in 1:9) {
    if (z[subco[i],j] < 0.01) {
      z[subco[i],j] <- 0.005
    }
    z[subco[i],j] <- log(z[subco[i],j]) - log(unsat[subco[i]])
  }
}

sdz <- apply(z[sam[which(sam %in% subco)],],2,sd)
z <- z %*% diag(1/sdz) # normalize by subcohort sd's

x <- data.matrix(mydat[,c(20, 23, 24, 25, 26)]) # auxiliary variables: dietary intake
x <- log(1+x) # log-transformation
x <- cbind(rep(1,nrow(x)), x) # add intercept term

Ysub <- Y[subco]
deltasub <- delta[subco]
zsub <- z[subco,]
N <- length(subco)

# at-risk indicators
risk <- matrix(0, nrow = length(fullcases), ncol = n)
for (i in 1:length(fullcases)) {
  for (j in 1:n) {
    if (Y[j] >= Y[fullcases[i]]) {
      risk[i,j] <- 1
    }
  }
}

# posterior predictive parameters
xihat <-solve(crossprod(cbind(x[subco,], w[subco,]))) %*% crossprod(cbind(x[subco,], w[subco,]), z[subco,])
Psi <- crossprod(z[subco,]-cbind(x[subco,], w[subco,]) %*% xihat)
VTVinv <- solve(crossprod(cbind(x[subco,], w[subco,]))) # inverse of V_S cross product
Ve <- eigen(VTVinv, symmetric = TRUE)
Vroot <- Ve$vectors %*% diag(sqrt(Ve$values)) %*% t(Ve$vectors) # square root of VTVinv, for generating matrix normal

# set starting beta values and proposal variance
current <- init_est
propvar <- (2.4^2)*(1/(dim(z)[2]+dim(w)[2]))* V_est

storesub <- exp(zsub %*% current[1:9]) # pre-calculate exponential terms for individuals with full covariates
Sigma <- riwish(N, Psi) # generate initial value of Sigma
numrow <- dim(w)[2]+dim(x)[2]

e <- eigen(Sigma, symmetric = TRUE)
root <- e$vectors %*% diag(sqrt(e$values)) %*% t(e$vectors) # square root of Sigma, for generating matrix normal

U_xi_curr <- matrix(rnorm(numrow*9), nrow = numrow, ncol = 9)
xi_curr <- xihat + Vroot %*% U_xi_curr %*% root # generate initial value of xi from matrix normal distribution

U_curr <- mvrnorm(num_miss, mu = rep(0,9), Sigma = diag(9))

Zmis <- U_curr %*% root + cbind(x[-subco,], w[-subco,]) %*% xi_curr
store <- numeric(n)
store[-subco] <- exp(Zmis %*% current[1:9])*exp(w[-subco,]  %*%  current[10:16])
store[subco] <- storesub*exp(w[subco,]  %*%  current[10:16])
partprop <- matrix(, nrow = n, ncol =1)
partprop[,1] <- store

denom <- risk %*% partprop # denominator factors of cox partial likelihood

hcurr <- prod(partprop[fullcases,1]*(n/2)/denom[,1]) # compute likelihood (n/2 factors due to numerical evaluation)

post_samp <- matrix(, nrow = iter, ncol = (dim(z)[2]+dim(w)[2])) # posterior samples of beta

# computes t-prior density value for a value of beta (3 degrees of freedom)
tprior <- function(beta) {
  return(prod(dt(beta,df=3)))
}

# Metropolis-Hastings
for (t in 1:iter) {
  
  prop <- mvrnorm(1, mu = current, Sigma = propvar) # generate beta proposal
  
  storesub <- exp(zsub %*% prop[1:9])
  
  e <- eigen(Sigma, symmetric = TRUE) # generate Sigma proposal
  root <- e$vectors %*% diag(sqrt(e$values)) %*% t(e$vectors) 
  U_xi_prop <- rhoxi* U_xi_curr + sqrt(1-rhoxi^2) *matrix(rnorm(numrow*9), nrow = numrow, ncol = 9) # generate xi proposal correlated with current value
  xi_prop <- xihat + Vroot %*% U_xi_prop %*% root

  U_prop <- rho*U_curr + sqrt(1-rho^2) *mvrnorm(num_miss, mu = rep(0,9), Sigma = diag(9)) # generate U proposal correlated with current value
  Zmis <- U_prop %*% root + cbind(x[-subco,], w[-subco,]) %*% xi_prop
  store <- numeric(n)
  store[-subco] <- exp(Zmis %*% prop[1:9])*exp(w[-subco,]  %*%  prop[10:16])
  store[subco] <- storesub*exp(w[subco,]  %*%  prop[10:16])
  partprop[,1] <- store
  
  
  denom <- risk %*% partprop
  
  hprop <- prod(partprop[fullcases,1]*(n/2)/denom[,1]) # compute proposal likelihood value
  
  if (hprop > 0) {
    ratio <- (hprop/hcurr)*(tprior(prop)/tprior(current)) # acceptance ratio
    a <- rbinom(n=1, 1, min(ratio,1))
    if (a == 1) {
      current <- prop
      hcurr <- hprop
      U_curr <- U_prop
      U_xi_curr <- U_xi_prop
    }
  }
  
  post_samp[t,] <- current
}

saveRDS(post_samp, file = "app_samp.rds")
end <- Sys.time()-start_time
print(end)




