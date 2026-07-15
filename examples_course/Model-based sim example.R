#########################################################################
# Advanced Small Area Estimation 
#
# Author:    Prof. Dr. Timo Schmid
# Content:   Model-based simulation - Point estimation                   
#########################################################################
rm(list=ls())

### Load required packages
library(emdi)
library(saeSim)
library(ggplot2)

#########################################################################
### Step 1: Data generation with R package -------------------------------------------------
#########################################################################

### a) Define input parameters and setup -------------------------------------------------
set.seed(100)
simruns <- 100  # Number of Monte-Carlo replications
Domains <- 50  # Number of Small Areas
# Define area-specific sample sizes
sample_size <- c(12, 23, 28, 14, 10, 23, 19, 25, 29, 10, 14, 18, 15, 20, 13, 12, 16, 
                 27, 20, 26, 27, 23, 12, 12, 11, 18, 17, 29, 11, 29, 17, 9, 14, 8, 
                  8, 18, 21, 21, 16, 16, 25, 13, 26, 19, 28, 20, 24, 9, 25, 21)
# Define area-specific population sizes
pop_size <- rep(200,Domains)
# First summary statistics
sum(pop_size) # Size of the population
sum(sample_size) # Size of the sample
summary(sample_size) # Summary of the area-specific sample sizes

# Help function (in order to keep the population in the data generation process)
calc_keepPop <- function(dat)
{
  attr(dat, "pop") = dat
  dat
}

### b) Define the data generation mechanism -------------------------------------------------
# Normal setup compare Rojas-Perilla et al. (2020) Table 6 
# Help functions for x,e and sampling
gen_XNorm <- function(dat, m = dat$muD, s = 3) {
  dat["x"] = rnorm(nrow(dat), mean = m, sd = s)
  return(dat)
}
gen_myE <- function(dat, m = 0, s = 1000) {
  dat["e"] = rnorm(nrow(dat), mean = m, sd = s)
  dat
}
sampler <- function(DAT) {
  smp = as.data.frame(matrix(nrow = sum(sample_size) , ncol = ncol(DAT)))
  brd = append(0, cumsum(sample_size))
  for (i in 1:Domains) {
    smp[((brd[i] + 1):brd[i + 1]), ] = (DAT[DAT$idD == i, ])[sample(1:sum(DAT$idD ==
                                                                            i), size = sample_size[i]), ]
  }
  attr(smp, "pop") = DAT
  colnames(smp) = colnames(DAT)
  return(smp)
}

# Define the components of one arbitrary population
setup <- sim_base(data = base_id(nDomains=Domains, nUnits=pop_size)) %>%
  sim_gen(gen_generic(runif, min = -3, max = 3, groupVars="idD", name = "muD")) %>%
  sim_gen(gen_XNorm)       %>%  
  as.data.frame %>%
  sim_gen(generator=gen_myE)         %>% 
  sim_gen_v(mean=0,sd=500)         %>% 
  sim_resp_eq(y = 4500 - 400*x + v + e)  

# Generate "simruns" populations defined by the setup
Pop <- sim(setup,R=simruns)
# Set the negative values to 0 (to avoid negative income data)
Pop <- lapply(Pop, function(X){ X$y[X$y<0] = 0; return (X) })
# Generate "simruns" samples from the population with "sample_size" area-specific sample sizes
Pop <- lapply(Pop,sampler)

#########################################################################
### Step 2: Computations on the true population --------------------------------------------
#########################################################################

# Define matrices for the true values of the indicators of interest
true_mean<-matrix(NA,Domains,simruns)
true_hcr<-matrix(NA,Domains,simruns)
true_pgap<-matrix(NA,Domains,simruns)
true_qsr<-matrix(NA,Domains,simruns)

# Help functions
hcr_function <- function(y, threshold) {
  mean(y < threshold)
}
qsr_function <- function(y) {
  sum(y[(y > quantile(y, 0.8))]) / sum(y[(y < quantile(y, 0.2))])
}
pgap_function <- function(y, threshold) {
  mean((y < threshold) * (threshold - y) / threshold)
}

# Calculate true values
for(i in 1:simruns){
  test<-attr(Pop[[i]],"pop")
  threshold<-0.6*median(test$y)
  
  true_mean[,i]<-tapply(test$y,test$idD,mean)
  true_hcr[,i]<-tapply(test$y,test$idD,hcr_function,threshold=threshold)
  true_qsr[,i]<-tapply(test$y,test$idD,qsr_function)
  true_pgap[,i]<-tapply(test$y,test$idD,pgap_function,threshold=threshold)
}

#########################################################################
### Step 3: Computations on the sample -----------------------------------------------------
#########################################################################

### a) Fit the models of interest -------------------------------------------------

# Fit the EBP without transformation 
EBP_raw <- NULL
for (i in 1:simruns) {
  test <- attr(Pop[[i]], "pop")
  EBP_raw[[i]] <- ebp(
    fixed = y ~ x,
    transformation = "no",
    pop_data = test,
    pop_domains = "idD",
    smp_data = Pop[[i]],
    smp_domains = "idD",
    MSE = FALSE,
    na.rm = TRUE
  )
  print(i)
}

# Fit the EBP with a fixed Log transformation 
EBP_log <- NULL
for (i in 1:simruns) {
  test <- attr(Pop[[i]], "pop")
  EBP_log[[i]] <- ebp(
    fixed = y ~ x,
    transformation = "log",
    pop_data = test,
    pop_domains = "idD",
    smp_data = Pop[[i]],
    smp_domains = "idD",
    MSE = FALSE,
    na.rm = TRUE
  )
  print(i)
}

# Fit the EBP with Box-Cox transformation 
EBP_boxcox <- NULL
for (i in 1:simruns) {
  test <- attr(Pop[[i]], "pop")
  EBP_boxcox[[i]] <- ebp(
    fixed = y ~ x,
    pop_data = test,
    pop_domains = "idD",
    smp_data = Pop[[i]],
    smp_domains = "idD",
    MSE = FALSE,
    na.rm = TRUE
  )
  print(i)
}

# Fit the EBP with a data-driven Log-shift transformation 
EBP_logshift <- NULL
for (i in 1:simruns) {
  test <- attr(Pop[[i]], "pop")
  EBP_logshift[[i]] <-
    ebp(
      fixed = y ~ x,
      transformation = "log.shift",
      pop_data = test,
      pop_domains = "idD",
      smp_data = Pop[[i]],
      smp_domains = "idD",
      MSE = FALSE,
      na.rm = TRUE
    )
  print(i)
}

# Fit the EBP with a dual-power transformation 
EBP_dual <- NULL
for (i in 1:simruns) {
  test <- attr(Pop[[i]], "pop")
  EBP_dual[[i]] <- ebp(
    fixed = y ~ x,
    transformation = "dual",
    pop_data = test,
    pop_domains = "idD",
    smp_data = Pop[[i]],
    smp_domains = "idD",
    MSE = FALSE,
    na.rm = TRUE
  )
  print(i)
}

### b) Extract the indicators of interest -------------------------------------------------
# Define matrices for the estimates of the indicators of interest
EBPest_raw_mean <- matrix(NA, Domains, simruns)
EBPest_raw_hcr <- matrix(NA, Domains, simruns)
EBPest_raw_pgap <- matrix(NA, Domains, simruns)
EBPest_raw_qsr <- matrix(NA, Domains, simruns)

EBPest_log_mean <- matrix(NA, Domains, simruns)
EBPest_log_hcr <- matrix(NA, Domains, simruns)
EBPest_log_pgap <- matrix(NA, Domains, simruns)
EBPest_log_qsr <- matrix(NA, Domains, simruns)

EBPest_boxcox_mean <- matrix(NA, Domains, simruns)
EBPest_boxcox_hcr <- matrix(NA, Domains, simruns)
EBPest_boxcox_pgap <- matrix(NA, Domains, simruns)
EBPest_boxcox_qsr <- matrix(NA, Domains, simruns)

EBPest_logshift_mean <- matrix(NA, Domains, simruns)
EBPest_logshift_hcr <- matrix(NA, Domains, simruns)
EBPest_logshift_pgap <- matrix(NA, Domains, simruns)
EBPest_logshift_qsr <- matrix(NA, Domains, simruns)

EBPest_dual_mean <- matrix(NA, Domains, simruns)
EBPest_dual_hcr <- matrix(NA, Domains, simruns)
EBPest_dual_pgap <- matrix(NA, Domains, simruns)
EBPest_dual_qsr <- matrix(NA, Domains, simruns)

# Fill the matrices with estimates
for(i in 1:simruns) {
  EBPest_raw_mean[, i] <- EBP_raw[[i]]$ind$Mean
  EBPest_raw_hcr[, i] <- EBP_raw[[i]]$ind$Head_Count
  EBPest_raw_pgap[, i] <- EBP_raw[[i]]$ind$Poverty_Gap
  EBPest_raw_qsr[, i] <- EBP_raw[[i]]$ind$Quintile_Share
  
  EBPest_log_mean[, i] <- EBP_log[[i]]$ind$Mean
  EBPest_log_hcr[, i] <- EBP_log[[i]]$ind$Head_Count
  EBPest_log_pgap[, i] <- EBP_log[[i]]$ind$Poverty_Gap
  EBPest_log_qsr[, i] <- EBP_log[[i]]$ind$Quintile_Share
  
  EBPest_boxcox_mean[, i] <- EBP_boxcox[[i]]$ind$Mean
  EBPest_boxcox_hcr[, i] <- EBP_boxcox[[i]]$ind$Head_Count
  EBPest_boxcox_pgap[, i] <- EBP_boxcox[[i]]$ind$Poverty_Gap
  EBPest_boxcox_qsr[, i] <- EBP_boxcox[[i]]$ind$Quintile_Share
  
  EBPest_logshift_mean[, i] <- EBP_logshift[[i]]$ind$Mean
  EBPest_logshift_hcr[, i] <- EBP_logshift[[i]]$ind$Head_Count
  EBPest_logshift_pgap[, i] <- EBP_logshift[[i]]$ind$Poverty_Gap
  EBPest_logshift_qsr[, i] <- EBP_logshift[[i]]$ind$Quintile_Share
  
  EBPest_dual_mean[, i] <- EBP_dual[[i]]$ind$Mean
  EBPest_dual_hcr[, i] <- EBP_dual[[i]]$ind$Head_Count
  EBPest_dual_pgap[, i] <- EBP_dual[[i]]$ind$Poverty_Gap
  EBPest_dual_qsr[, i] <- EBP_dual[[i]]$ind$Quintile_Share
}

#########################################################################
### Step 4: Assess the performance by quality measures -------------------------------------
#########################################################################

# Help function with quality measures
Qualitymeasures <- function (True.ind, Est.ind) {
  True.mean <- True.ind
  Est.mean <- Est.ind
  m <- dim(True.mean)[1]
  NoSim <- dim(True.mean)[2]
  True.mean <- t(True.mean)
  Est.mean <- t(Est.mean)
  
  
  assign("m", m, pos = 1)
  Bias <- rep(0, m)
  RB <- rep(0, m)
  RRMSE <- rep(0, m)
  
  True.MSE <- rep(0, m)
  Est.MSE <- rep(0, m)
  CV.MSE <- rep(0, m)
  BIAS.MSE <- rep(0, m)
  
  for (i in 1:m) {
    RB[i] <- mean(((Est.mean[, i]) - (True.mean[, i])) / (True.mean[, i]))
    Bias[i] <- mean((Est.mean[, i]) - (True.mean[, i]))
    RRMSE[i] <-
      sqrt(mean((Est.mean[, i] - True.mean[, i]) ^ 2)) / mean(True.mean[, i])
    True.MSE[i] <- mean((Est.mean[, i] - True.mean[, i]) ^ 2)
  }
  True.RootMSE = sqrt(True.MSE)
  
  list(
    Bias = Bias,
    RelBias = RB,
    RelRMSE = RRMSE,
    EmpRMSE = True.RootMSE
  )
}

# Example: Calculate the quality measures for EBP Raw for indicator mean
Qualitymeasures(True.ind = true_mean, Est.ind = EBPest_raw_mean)

#########################################################################
### Step 5: Present the quality measures ---------------------------------------------------
#########################################################################
# Extract results for QSR
EBP_raw_rmse <-
  Qualitymeasures(True.ind = true_qsr, Est.ind = EBPest_raw_qsr)$EmpRMSE
EBP_raw_bias <-
  Qualitymeasures(True.ind = true_qsr, Est.ind = EBPest_raw_qsr)$Bias
EBP_log_rmse <-
  Qualitymeasures(True.ind = true_qsr, Est.ind = EBPest_log_qsr)$EmpRMSE
EBP_log_bias <-
  Qualitymeasures(True.ind = true_qsr, Est.ind = EBPest_log_qsr)$Bias
EBP_boxcox_rmse <-
  Qualitymeasures(True.ind = true_qsr, Est.ind = EBPest_boxcox_qsr)$EmpRMSE
EBP_boxcox_bias <-
  Qualitymeasures(True.ind = true_qsr, Est.ind = EBPest_boxcox_qsr)$Bias
EBP_logshift_rmse <-
  Qualitymeasures(True.ind = true_qsr, Est.ind = EBPest_logshift_qsr)$EmpRMSE
EBP_logshift_bias <-
  Qualitymeasures(True.ind = true_qsr, Est.ind = EBPest_logshift_qsr)$Bias
EBP_dual_rmse <-
  Qualitymeasures(True.ind = true_qsr, Est.ind = EBPest_dual_qsr)$EmpRMSE
EBP_dual_bias <-
  Qualitymeasures(True.ind = true_qsr, Est.ind = EBPest_dual_qsr)$Bias

# Bias
BiasEbp <-
  as.data.frame(c(
    EBP_raw_bias,
    EBP_log_bias,
    EBP_boxcox_bias,
    EBP_logshift_bias,
    EBP_dual_bias
  ))
names(BiasEbp)[names(BiasEbp) == "c(EBP_raw_bias, EBP_log_bias, EBP_boxcox_bias, EBP_logshift_bias, EBP_dual_bias)"] <-
  "bias"
BiasEbp$method <- 0
BiasEbp$method[1:Domains] <- "EBP"
BiasEbp$method[(Domains + 1):(2 * Domains)] <- "EBP Log"
BiasEbp$method[(2 * Domains + 1):(3 * Domains)] <- "EBP Box Cox"
BiasEbp$method[(3 * Domains + 1):(4 * Domains)] <- "EBP Log-Shift"
BiasEbp$method[(4 * Domains + 1):(5 * Domains)] <- "EBP Dual"
ggplot(BiasEbp, aes(x = method, y = bias, fill = method)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, color = "red") +
  labs(y = "Bias", x = "", title = "QSR") +
  theme(legend.position = "none")

# RMSE
RmseEbp <-
  as.data.frame(c(
    EBP_raw_rmse,
    EBP_log_rmse,
    EBP_boxcox_rmse,
    EBP_logshift_rmse,
    EBP_dual_rmse
  ))
names(RmseEbp)[names(RmseEbp) == "c(EBP_raw_rmse, EBP_log_rmse, EBP_boxcox_rmse, EBP_logshift_rmse, EBP_dual_rmse)"] <-
  "rmse"
RmseEbp$method <- 0
RmseEbp$method[1:Domains] <- "EBP"
RmseEbp$method[(Domains + 1):(2 * Domains)] <- "EBP Log"
RmseEbp$method[(2 * Domains + 1):(3 * Domains)] <- "EBP Box Cox"
RmseEbp$method[(3 * Domains + 1):(4 * Domains)] <- "EBP Log-Shift"
RmseEbp$method[(4 * Domains + 1):(5 * Domains)] <- "EBP Dual"
ggplot(RmseEbp, aes(x = method, y = rmse, fill = method)) +
  geom_boxplot() +
  labs(y = "RMSE", x = "", title = "QSR") +
  theme(legend.position = "none")

