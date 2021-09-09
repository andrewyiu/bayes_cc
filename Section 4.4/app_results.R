##############################################################################################
## Create the density strips plot in Figure 2 and the "Analysis results" section of Table 2 ##
##############################################################################################

library(denstrip)

# load application results
results <- readRDS("app_samp.rds")

# remove warm-up iteraions and confounder parameters
results <- results[200001:1000000,1:9]

col1 <- "#FE9929"; col2="#4292C6"; col3="#41AB5D"; col4="#6A51A3"
colmin <- "gray92"
wd <- 0.3
ci <- c(0.025, 0.975)
base <- 1
dy <- 0.5
dg <- 0.7
base2 <- base + 4*dy + dg
base3 <- base2 + 3*dy + dg

prev.plot <- function(sam) { 
  
  par(mar=c(3, 0, 0, 0), mgp=c(2,1,0))
  
  xmax <- 1.8
  plot(0, type="n", xlim=c(0.4, xmax), ylim=c(1,7), axes=FALSE, xlab="Hazard ratio", ylab="")
  lim <- par("usr")
  rect(0.6, lim[3], 1.7, lim[4], col="gray92", border="gray92")
  axis(1, at=seq(0.6,1.7,by=0.1))
  abline(v=seq(0.6,1.7,by=0.1), col="white")
  
  denstrip(sam[,1],        at=base3+dy,  width=wd, colmax=col1, colmin=colmin, from=0.6, to=1.7, ticks=quantile(sam[,1], ci))
  denstrip(sam[,2],   at=base3, width=wd, colmax=col2, colmin=colmin, from=0.6, to=1.7, ticks=quantile(sam[,2], ci))
  
  denstrip(sam[,3], at=base2+2*dy, width=wd, colmax=col1, colmin=colmin, from=0.6, to=1.7, ticks=quantile(sam[,3], ci))
  denstrip(sam[,4],        at=base2+dy,  width=wd, colmax=col2, colmin=colmin, from=0.6, to=1.7, ticks=quantile(sam[,4], ci))
  denstrip(sam[,5],   at=base2, width=wd, colmax=col3, colmin=colmin, from=0.6, to=1.7, ticks=quantile(sam[,5], ci))
  
  denstrip(sam[,6], at=base+3*dy, width=wd, colmax=col1, colmin=colmin, from=0.6, to=1.7, ticks=quantile(sam[,6], ci))
  denstrip(sam[,7],        at=base+2*dy, width=wd, colmax=col2, colmin=colmin, from=0.6, to=1.7, ticks=quantile(sam[,7], ci))
  denstrip(sam[,8],   at=base+dy,  width=wd, colmax=col3, colmin=colmin, from=0.6, to=1.7, ticks=quantile(sam[,8], ci))
  denstrip(sam[,9], at=base,  width=wd, colmax=col4, colmin=colmin, from=0.6, to=1.7, ticks=quantile(sam[,9], ci))
  
  text(0.350, c((base3+2*dy), (base2+3*dy)), c("ocSFAs", "ecSFAs"), pos=4)
  text(0.350, c(base + 4*dy), c("vlcSFAs"), pos=4)
  # labs <- c("Overall","Diagnosed","Undiagnosed")
  # labsn <- c(expression(paste("Overall ",pi[N])), expression(paste("Diagnosed ", (pi*delta)[N])), expression(paste("Undiagnosed ", bar((pi*delta))[N])))
  # labsg <- c(expression(paste("Overall ",pi[G])), expression(paste("Diagnosed ", (pi*delta)[G])), expression(paste("Undiagnosed ", bar((pi*delta))[G])))
  text(0.400, base+c(3*dy,2*dy,dy,0), c("C20:0", "C22:0", "C23:0", "C24:0"), pos=4, col=c(col1,col2,col3, col4))
  text(0.40, base2+c(2*dy,dy,0), c("C14:0", "C16:0", "C18:0"), pos=4, col=c(col1,col2,col3))
  text(0.40, base3+c(dy,0), c("C15:0", "C17:0"), pos=4, col=c(col1,col2))
  
}

## create Figure 4
prev.plot(exp(results))
abline(v = 1.00, col = "black", lty = 5)

## create "Analysis results" section of Table 4
# take exponential transformation
results <- exp(results)

hr_mean <- round(colMeans(results),2) # posterior means
int_low <- round(apply(results, 2,FUN=quantile, prob = 0.025),2) # lower 2.5% credible quantiles
int_high <- round(apply(results, 2,FUN=quantile, prob = 0.975),2) # upper 2.5% credible quantiles
# posterior probabilities of HR <= 1
pr1 <- numeric(9)
for (i in 1:9) {
  pr1[i] <- round(length(which(results[,i] <= 1))/800000,3)
}

SFA <- c("C15:0", "C17:0", "C14:0", "C16:0", "C18:0", "C20:0", "C22:0", "C23:0", "C24:0")
tab <- data.frame(cbind(SFA, hr_mean, int_low, int_high, pr1))
print(tab)
