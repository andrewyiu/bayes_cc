#################################
## Create Figure 2 and Table 3 ##
#################################

library(survival)
library(ggplot2)
library(gridExtra)
library(reshape)

# load data
data = read.csv("mydata_up.csv")
for (j in 1:25653) {
  if (is.na(data$dmstatus_ver_outc[j]) == TRUE) {
    data$dmstatus_ver_outc[j] <- 3
  }
}

# load synthetic data experiment results
results <- readRDS("synth_results.rds")

# remove prevalent cases
mydat_orig  <- data[-which(data$dmstatus_ver == 2),]

# remove individuals with unknown status
mydat_orig  <- mydat_orig[-which(mydat_orig$dmstatus_ver == -1),]

# remove individuals with missing confounders
mydat_orig  <- mydat_orig[which(is.na(mydat_orig$pa_index)==FALSE),]
mydat_orig  <- mydat_orig[which(is.na(mydat_orig$bmi_adj)==FALSE),]
mydat_orig  <- mydat_orig[which(is.na(mydat_orig$waist_adj)==FALSE),]
mydat_orig  <- mydat_orig[which(is.na(mydat_orig$sex)==FALSE),]
mydat_orig  <- mydat_orig[which(is.na(mydat_orig$age_recr)==FALSE),]

# remove excleier variables
mydat_orig  <- mydat_orig[-which(mydat_orig$excleier == 1),]
mydat_orig  <- mydat_orig[-which(mydat_orig$excleier == 3),]
mydat_orig  <- mydat_orig[-which(is.na(mydat_orig$excleier)==TRUE),]

#remove cases and subcohort members who do not have fatty acid measurements
mydat_orig  <- mydat_orig[-which(mydat_orig$outcome==1 & is.na(mydat_orig$C100)==TRUE),]
mydat_orig  <- mydat_orig[-which(mydat_orig$outcome==0 & mydat_orig$dmstatus_ver == 1 & is.na(mydat_orig$C100)==TRUE),]

subcohort <- which(mydat_orig$dmstatus_ver_outc<=1) # original subcohort

set.seed(123)
n_big <- 2000000
impute <- sample(subcohort, n_big, replace = TRUE) 
mydat <- mydat_orig[impute,]
Y <- mydat$fup_time
delta <- as.integer(mydat$dmstatus_ver_outc == 1 | mydat$dmstatus_ver_outc == 2)

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
w <- w %*% diag((1/c(1,sdw[2:4],1,1,1)))

# additive logratio transformation (replace any zero fatty acid measurements with 0.005 first)
for (i in 1:n_big) {
  for (j in 1:9) {
    if (z[i,j] < 0.01) {
      
      z[i,j] <- 0.005
      
    }
    z[i,j] <- log(z[i,j]) - log(unsat[i])
    
  }
}

z <- z %*% diag((1/apply(z,2,sd))) # normalize by sd's


testfull <- list(time= as.numeric(Y),
                 status=delta, C15 = z[,1],
                 C17 = z[,2],
                 C14 = z[,3],
                 C16 = z[,4],
                 C18 = z[,5],
                 C20 = z[,6],
                 C22 = z[,7],
                 C23 = z[,8],
                 C24 = z[,9],
                 sex = w[,1],
                 age = w[,2],
                 waist = w[,3],
                 bmi = w[,4],
                 edu1 = w[,5],
                 edu2 = w[,6],
                 edu3 = w[,7]
)
fitfull <- coxph(Surv(time,status) ~ C15 +
                   C17 +
                   C14 +
                   C16 +
                   C18 +
                   C20 +
                   C22 +
                   C23 +
                   C24 +
                   sex +
                   age +
                   waist +
                   bmi +
                   edu1 +
                   edu2 +
                   edu3
                 ,data=testfull, ties = "breslow")

truth <-fitfull$coefficients # monte carlo estimate of "true" beta values

## plot violin plots in figure 2
i <- 1

av <- cbind(results[,i],results[,(16+i)])
ave <- as.vector(av)
grps <- c(rep("Bayes", dim(results)[1]),rep("Prentice", dim(results)[1]))
d <- cbind(ave, grps)
d <- data.frame(d)
df <- melt(d, id.vars = "grps")

plot1 <- ggplot(df, aes(x=factor(df$grps), y=ave, fill = factor(df$grps))) + 
  geom_violin() + scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  theme_light() + theme(axis.title.x=element_blank(),
                        axis.ticks.x=element_blank()) +
  ylab("log-hazard ratio")  + geom_hline(yintercept=truth[i], linetype="dashed", color = "red", size = 1.0) + 
  ggtitle("C15:0")+ theme(legend.position='none')

i <- i + 1

av <- cbind(results[,i],results[,(16+i)])
ave2 <- as.vector(av)
grps <- c(rep("Bayes", dim(results)[1]),rep("Prentice", dim(results)[1]))
d <- cbind(ave, grps)
d<- data.frame(d)
df <- melt(d, id.vars = "grps")

plot2 <- ggplot(df, aes(x=factor(df$grps), y=ave2, fill = factor(df$grps))) + 
  geom_violin() + scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  theme_light() + theme(axis.title.x=element_blank(),
                        axis.ticks.x=element_blank()) +
  ylab("log-hazard ratio")  + geom_hline(yintercept=truth[i], linetype="dashed", color = "red", size = 1) + 
  ggtitle("C17:0")+ theme(legend.position='none')

i <- i + 1

av <- cbind(results[,i],results[,(16+i)])
ave3 <- as.vector(av)
grps <- c(rep("Bayes", dim(results)[1]),rep("Prentice", dim(results)[1]))
d <- cbind(ave, grps)
d <- data.frame(d)
df <- melt(d, id.vars = "grps")

plot3 <- ggplot(df, aes(x=factor(grps), y=ave3, fill = factor(grps))) + 
  geom_violin() + scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  theme_light() + theme(axis.title.x=element_blank(),
                        axis.ticks.x=element_blank()) +
  ylab("log-hazard ratio")  + geom_hline(yintercept=truth[i], linetype="dashed", color = "red", size = 1) + 
  ggtitle("C14:0")+ theme(legend.position='none')

i <- i + 1

av <- cbind(results[,i],results[,(16+i)])
ave4 <- as.vector(av)
grps <- c(rep("Bayes", dim(results)[1]),rep("Prentice", dim(results)[1]))
d <- cbind(ave, grps)
d <- data.frame(d)
df <- melt(d, id.vars = "grps")

plot4 <- ggplot(df, aes(x=factor(grps), y=ave4, fill = factor(grps))) + 
  geom_violin() + scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  theme_light() + theme(axis.title.x=element_blank(),
                        axis.ticks.x=element_blank()) +
  ylab("log-hazard ratio")  + geom_hline(yintercept=truth[i], linetype="dashed", color = "red", size = 1) + 
  ggtitle("C16:0")+ theme(legend.position='none')
i <- i + 1

av <- cbind(results[,i],results[,(16+i)])
ave5 <- as.vector(av)
grps <- c(rep("Bayes", dim(results)[1]),rep("Prentice", dim(results)[1]))
d <- cbind(ave, grps)
d <- data.frame(d)
df <- melt(d, id.vars = "grps")

plot5 <- ggplot(df, aes(x=factor(grps), y=ave5, fill = factor(grps))) + 
  geom_violin() + scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  theme_light() + theme(axis.title.x=element_blank(),
                        axis.ticks.x=element_blank()) +
  ylab("log-hazard ratio")  + geom_hline(yintercept=truth[i], linetype="dashed", color = "red", size = 1) + 
  ggtitle("C18:0")+ theme(legend.position='none')
i <- i + 1

av <- cbind(results[,i],results[,(16+i)])
ave6 <- as.vector(av)
grps <- c(rep("Bayes", dim(results)[1]),rep("Prentice", dim(results)[1]))
d <- cbind(ave, grps)
d <- data.frame(d)
df <- melt(d, id.vars = "grps")

plot6 <- ggplot(df, aes(x=factor(grps), y=ave6, fill = factor(grps))) + 
  geom_violin() + scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  theme_light() + theme(axis.title.x=element_blank(),
                        axis.ticks.x=element_blank()) +
  ylab("log-hazard ratio")  + geom_hline(yintercept=truth[i], linetype="dashed", color = "red", size = 1) + 
  ggtitle("C20:0")+ theme(legend.position='none')
i <- i + 1

av <- cbind(results[,i],results[,(16+i)])
ave7 <- as.vector(av)
grps <- c(rep("Bayes", dim(results)[1]),rep("Prentice", dim(results)[1]))
d <- cbind(ave, grps)
d <- data.frame(d)
df <- melt(d, id.vars = "grps")

plot7 <- ggplot(df, aes(x=factor(grps), y=ave7, fill = factor(grps))) + 
  geom_violin() + scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  theme_light() + theme(axis.title.x=element_blank(),
                        axis.ticks.x=element_blank()) +
  ylab("log-hazard ratio")  + geom_hline(yintercept=truth[i], linetype="dashed", color = "red", size = 1) + 
  ggtitle("C22:0")+ theme(legend.position='none') 
i <- i + 1

av <- cbind(results[,i],results[,(16+i)])
ave8 <- as.vector(av)
grps <- c(rep("Bayes", dim(results)[1]),rep("Prentice", dim(results)[1]))
d <- cbind(ave, grps)
d <- data.frame(d)
df <- melt(d, id.vars = "grps")

plot8 <- ggplot(df, aes(x=factor(grps), y=ave8, fill = factor(grps))) + 
  geom_violin() + scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  theme_light() + theme(axis.title.x=element_blank(),
                        axis.ticks.x=element_blank()) +
  ylab("log-hazard ratio")  + geom_hline(yintercept=truth[i], linetype="dashed", color = "red", size = 1) + 
  ggtitle("C23:0")+ theme(legend.position='none')
i <- i + 1

av <- cbind(results[,i],results[,(16+i)])
ave9 <- as.vector(av)
grps <- c(rep("Bayes", dim(results)[1]),rep("Prentice", dim(results)[1]))
d <- cbind(ave, grps)
d <- data.frame(d)
df <- melt(d, id.vars = "grps")

plot9 <- ggplot(df, aes(x=factor(grps), y=ave9, fill = factor(grps))) + 
  geom_violin() + scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  theme_light() + theme(axis.title.x=element_blank(),
                        axis.ticks.x=element_blank()) +
  ylab("log-hazard ratio")  + geom_hline(yintercept=truth[i], linetype="dashed", color = "red", size = 1) + 
  ggtitle("C24:0")+ theme(legend.position='none')

grid.arrange(plot3, plot1, plot4,plot2,  plot5, plot6, plot7, plot8, plot9, nrow = 3, ncol = 3)
## saving the plot
#g <- arrangeGrob(plot3, plot1, plot4,plot2,  plot5, plot6, plot7, plot8, plot9, nrow = 3, ncol = 3)
#ggsave(file = "violin.pdf", g)

## computing the result in table 2
res <- results[,c(3,1,4,2,5,6:9,19,17,20,18,21,22:25)] # trim confounder variables and reorder

bias <- colMeans(res) - truth[c(3,1,4,2,5,6:9)]
stan_dev <- apply(res, 2, sd)
mse <- bias^2+stan_dev^2
rmse <- sqrt(mse)

bias <- round(c(rbind(bias[10:18],bias[1:9] )), 3)
stan_dev <- round(c(rbind(stan_dev[10:18],stan_dev[1:9])),3)
rmse <- round(c(rbind(rmse[10:18],rmse[1:9])),3)
group <- rep(c("Prentice", "Bayes"), 9)
SFA <- c(rep("C14", 2),rep("C15", 2),rep("C16", 2),rep("C17", 2),rep("C18", 2),rep("C20", 2),rep("C22", 2),rep("C23", 2),rep("C24", 2))
tab <- data.frame(cbind(group, SFA,bias, stan_dev,rmse)) 

## table 3
print(tab)
eg <- round(100*(mse[10:18]/mse[1:9])-100,0) ## efficiency gain
print(eg)
