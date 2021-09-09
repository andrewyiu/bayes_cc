########################################################################################
## Create the correlation heat map (Figure 1) and the data summary section of Table 2 ##
########################################################################################

library(ggplot2)
library(reshape)

# load data
data = read.csv("mydata_up.csv")

for (j in 1:25653) {
  if (is.na(data$dmstatus_ver_outc[j]) == TRUE) {
    data$dmstatus_ver_outc[j] <- 3
  }
}

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

comp <- which(is.na(mydat$C100)==FALSE) # individuals with no missing fatty acid measurements
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

cormat <- round(cor(z[comp[which(comp %in% sam)],]),2) # correlations before transformation
mean_raw <- round(colMeans(z[comp[which(comp %in% sam)],]),2) # means before transformation
sd_raw <- round(apply(z[comp[which(comp %in% sam)],], 2, sd),2) # sd's before transformation

# additive logratio transformation (replace any zero fatty acid measurements with 0.005 first)
for (i in 1:length(comp)) {
  for (j in 1:9) {
    if (z[comp[i],j] < 0.01) {
      z[comp[i],j] <- 0.005
    }
    z[comp[i],j] <- log(z[comp[i],j]) - log(unsat[comp[i]])
  }
}

cormat2 <- round(cor(z[comp[which(comp %in% sam)],]),2) # correlations after transformation
mean_alr <- round(colMeans(z[comp[which(comp %in% sam)],]),2) # means after transformation
sd_alr <- round(apply(z[comp[which(comp %in% sam)],], 2, sd),2) # sd's after transformation

cormat_full <- matrix(, nrow = 9, ncol = 9)
for (i in 1:9) {
  for (j in 1:9) {
    if (i < j) {
      cormat_full[i,j] <- cormat2[i,j]
    } else if (i > j) {
      cormat_full[i,j] <- cormat[i,j]
    } else {
      cormat_full[i,j] <- NA
    }
  }
}
colnames(cormat_full) <- c("C15:0", "C17:0", "C14:0", "C16:0", "C18:0", "C20:0", "C22:0", "C23:0", "C24:0")
rownames(cormat_full) <- c("C15:0", "C17:0", "C14:0", "C16:0", "C18:0", "C20:0", "C22:0", "C23:0", "C24:0")

melted_cormat_full <- melt(cormat_full)
melted_cormat_full$X1 <- factor(melted_cormat_full$X1, levels = c("C15:0", "C17:0", "C14:0", "C16:0", "C18:0", "C20:0", "C22:0", "C23:0", "C24:0"))
melted_cormat_full$X2 <- factor(melted_cormat_full$X2, levels = c("C15:0", "C17:0", "C14:0", "C16:0", "C18:0", "C20:0", "C22:0", "C23:0", "C24:0"))
ggplot(data = melted_cormat_full, aes(x=X1, y=X2, fill=value)) + 
  geom_tile() + scale_fill_gradient2(low = "#084594", high = "#CB181D", mid = "white", 
                                     midpoint = 0, limit = c(-1,1), space = "Lab", 
                                     name="Correlation") +
  theme_minimal() +
  geom_text(aes(label = sprintf("%0.2f", value)), data =  subset(melted_cormat_full,!is.na(value))) +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank()) 

SFA <- c("C15:0", "C17:0", "C14:0", "C16:0", "C18:0", "C20:0", "C22:0", "C23:0", "C24:0")
tab <- data.frame(cbind(SFA, mean_raw, sd_raw, mean_alr, sd_alr))
print(tab)