#########################################################################
# Advanced Small Area Estimation 
#
# Author:    Prof. Dr. Timo Schmid
# Content:   Model-based simulation - MSE estimation                   
#########################################################################
rm(list=ls())

### Load required packages
library(emdi)
library(saeSim)
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)

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

### a) Fit the EBP with Box-Cox transformation -------------------------------------------------

# Fit the EBP Box-Cox with a parametric bootstrap 
EBP_boxcox_para <- NULL
for (i in 1:simruns) {
  test <- attr(Pop[[i]], "pop")
  EBP_boxcox_para[[i]] <- ebp(
    fixed = y ~ x,
    pop_data = test,
    pop_domains = "idD",
    smp_data = Pop[[i]],
    smp_domains = "idD",
    boot_type = "parametric",
    B = 20, 
    MSE = TRUE,
    na.rm = TRUE
  )
  print(i)
}

# Fit the EBP Box-Cox with a wild bootstrap 
EBP_boxcox_wild <- NULL
for (i in 1:simruns) {
  test <- attr(Pop[[i]], "pop")
  EBP_boxcox_wild[[i]] <- ebp(
    fixed = y ~ x,
    pop_data = test,
    pop_domains = "idD",
    smp_data = Pop[[i]],
    smp_domains = "idD",
    boot_type = "wild",
    B = 20, 
    MSE = TRUE,
    na.rm = TRUE
  )
  print(i)
}



### b) Extract the indicators of interest -------------------------------------------------
# Define matrices for the estimates of the indicators of interest
EBPest_boxcox_mean <- matrix(NA, Domains, simruns)
EBPest_boxcox_hcr <- matrix(NA, Domains, simruns)
EBPest_boxcox_pgap <- matrix(NA, Domains, simruns)
EBPest_boxcox_qsr <- matrix(NA, Domains, simruns)

# Define matrices for the MSE estimates of the indicators of interest
EBPest_boxcox_para_mean <- matrix(NA, Domains, simruns)
EBPest_boxcox_para_hcr <- matrix(NA, Domains, simruns)
EBPest_boxcox_para_pgap <- matrix(NA, Domains, simruns)
EBPest_boxcox_para_qsr <- matrix(NA, Domains, simruns)

EBPest_boxcox_wild_mean <- matrix(NA, Domains, simruns)
EBPest_boxcox_wild_hcr <- matrix(NA, Domains, simruns)
EBPest_boxcox_wild_pgap <- matrix(NA, Domains, simruns)
EBPest_boxcox_wild_qsr <- matrix(NA, Domains, simruns)


# Fill the matrices with point and MSE estimates
for(i in 1:simruns) {
  EBPest_boxcox_mean[, i] <- EBP_boxcox_para[[i]]$ind$Mean
  EBPest_boxcox_hcr[, i] <- EBP_boxcox_para[[i]]$ind$Head_Count
  EBPest_boxcox_pgap[, i] <- EBP_boxcox_para[[i]]$ind$Poverty_Gap
  EBPest_boxcox_qsr[, i] <- EBP_boxcox_para[[i]]$ind$Quintile_Share
  
  EBPest_boxcox_para_mean[, i] <- EBP_boxcox_para[[i]]$MSE$Mean
  EBPest_boxcox_para_hcr[, i] <- EBP_boxcox_para[[i]]$MSE$Head_Count
  EBPest_boxcox_para_pgap[, i] <- EBP_boxcox_para[[i]]$MSE$Poverty_Gap
  EBPest_boxcox_para_qsr[, i] <- EBP_boxcox_para[[i]]$MSE$Quintile_Share
  
  EBPest_boxcox_wild_mean[, i] <- EBP_boxcox_wild[[i]]$MSE$Mean
  EBPest_boxcox_wild_hcr[, i] <- EBP_boxcox_wild[[i]]$MSE$Head_Count
  EBPest_boxcox_wild_pgap[, i] <- EBP_boxcox_wild[[i]]$MSE$Poverty_Gap
  EBPest_boxcox_wild_qsr[, i] <- EBP_boxcox_wild[[i]]$MSE$Quintile_Share
  
}

#########################################################################
### Step 4: Assess the performance by quality measures -------------------------------------
#########################################################################

# Help function with quality measures
Qualitymeasures_MSE <- function (True.ind, Est.ind, MSE) {
  True.mean <- True.ind
  Est.mean <- Est.ind
  m <- dim(True.mean)[1]
  NoSim <- dim(True.mean)[2]
  True.mean <- t(True.mean)
  Est.mean <- t(Est.mean)
  MSE <- t(MSE)
  
  
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
    Est.MSE[i]<-mean(MSE[,i])
    CV.MSE[i]<-(sqrt(mean((sqrt(MSE[,i])-sqrt(True.MSE[i]))^2))/sqrt(True.MSE[i]))*100
    BIAS.MSE[i]<-mean(MSE[,i]-True.MSE[i])/True.MSE[i]
  }
  True.RootMSE = sqrt(True.MSE)
  Est.RootMSE=sqrt(Est.MSE)
  RB.estRRMSE=((Est.RootMSE-True.RootMSE)/True.RootMSE)*100
  
  
  CR<-rep(0,m)
  CR_Length<-rep(0,m)
  coverage<-matrix(0,NoSim,m)
  Length<-matrix(0,NoSim,m)
  
  for (k in 1:NoSim)
  {
    for (oo in 1:m)
    {Length[k,oo]<-2*1.96*sqrt(MSE[k,oo])
    if ((Est.mean[k,oo]-1.96*sqrt(MSE[k,oo])< True.mean[k,oo] && True.mean[k,oo] < Est.mean[k,oo]+1.96*sqrt(MSE[k,oo]) ))coverage[k,oo]<-1
    else
      coverage[k,oo]<- 0}
  }
  
  for (i in 1:m){
    CR[i]<-mean(coverage[,i])
    CR_Length[i]<-mean(Length[,i])
  }
  
  list(
    Bias = Bias,
    RelBias = RB,
    RelRMSE = RRMSE,
    EmpRMSE = True.RootMSE,
    EstRMSE = Est.RootMSE,
    RelBias_RMSE = RB.estRRMSE,
    RelRMSE_RMSE = CV.MSE,
    Coverage = CR,
    CI_Length = CR_Length
  )
}

# Example: Calculate the quality measures for EBP Box-Cox with para. bootstrap for indicator mean
test<-Qualitymeasures_MSE(True.ind = true_mean, Est.ind = EBPest_boxcox_mean, MSE = EBPest_boxcox_para_mean)
summary(test$RelBias_RMSE)
        


#########################################################################
### Step 5: Present the quality measures ---------------------------------------------------
#########################################################################
# Extract results for the Mean
EBP_boxcox_mean_para_emprmse <-
  Qualitymeasures_MSE(True.ind = true_mean, Est.ind = EBPest_boxcox_mean, MSE = EBPest_boxcox_para_mean)$EmpRMSE
EBP_boxcox_mean_para_estrmse <-
  Qualitymeasures_MSE(True.ind = true_mean, Est.ind = EBPest_boxcox_mean, MSE = EBPest_boxcox_para_mean)$EstRMSE
EBP_boxcox_mean_wild_estrmse <-
  Qualitymeasures_MSE(True.ind = true_mean, Est.ind = EBPest_boxcox_mean, MSE = EBPest_boxcox_wild_mean)$EstRMSE

# Extract results for the HCR
EBP_boxcox_hcr_para_emprmse <-
  Qualitymeasures_MSE(True.ind = true_hcr, Est.ind = EBPest_boxcox_hcr, MSE = EBPest_boxcox_para_hcr)$EmpRMSE
EBP_boxcox_hcr_para_estrmse <-
  Qualitymeasures_MSE(True.ind = true_hcr, Est.ind = EBPest_boxcox_hcr, MSE = EBPest_boxcox_para_hcr)$EstRMSE
EBP_boxcox_hcr_wild_estrmse <-
  Qualitymeasures_MSE(True.ind = true_hcr, Est.ind = EBPest_boxcox_hcr, MSE = EBPest_boxcox_wild_hcr)$EstRMSE

# Extract results for the PGAP
EBP_boxcox_pgap_para_emprmse <-
  Qualitymeasures_MSE(True.ind = true_pgap, Est.ind = EBPest_boxcox_pgap, MSE = EBPest_boxcox_para_pgap)$EmpRMSE
EBP_boxcox_pgap_para_estrmse <-
  Qualitymeasures_MSE(True.ind = true_pgap, Est.ind = EBPest_boxcox_pgap, MSE = EBPest_boxcox_para_pgap)$EstRMSE
EBP_boxcox_pgap_wild_estrmse <-
  Qualitymeasures_MSE(True.ind = true_pgap, Est.ind = EBPest_boxcox_pgap, MSE = EBPest_boxcox_wild_pgap)$EstRMSE

# Extract results for the QSR
EBP_boxcox_qsr_para_emprmse <-
  Qualitymeasures_MSE(True.ind = true_qsr, Est.ind = EBPest_boxcox_qsr, MSE = EBPest_boxcox_para_qsr)$EmpRMSE
EBP_boxcox_qsr_para_estrmse <-
  Qualitymeasures_MSE(True.ind = true_qsr, Est.ind = EBPest_boxcox_qsr, MSE = EBPest_boxcox_para_qsr)$EstRMSE
EBP_boxcox_qsr_wild_estrmse <-
  Qualitymeasures_MSE(True.ind = true_qsr, Est.ind = EBPest_boxcox_qsr, MSE = EBPest_boxcox_wild_qsr)$EstRMSE



# Start to creat Figure in ggplot
# Domains
Domain <- 1:length(EBP_boxcox_mean_para_emprmse)

# ---- Mean Indicator ----
df_mean_long <- data.frame(
  Domain = Domain,
  EmpRMSE = EBP_boxcox_mean_para_emprmse,
  ParaEstRMSE = EBP_boxcox_mean_para_estrmse,
  WildEstRMSE = EBP_boxcox_mean_wild_estrmse
) %>%
  pivot_longer(cols = -Domain, names_to = "Method", values_to = "RMSE")

p_mean <- ggplot(df_mean_long, aes(x = Domain, y = RMSE, color = Method)) +
  geom_line(size = 1) +
  geom_point() +
  labs(title = "Mean Indicator", x = "Domain", y = "RMSE") +
  scale_color_manual(values = c("EmpRMSE" = "black", 
                                "ParaEstRMSE" = "red", 
                                "WildEstRMSE" = "orange")) +
  theme_gray()

# ---- HCR Indicator ----
df_hcr_long <- data.frame(
  Domain = Domain,
  EmpRMSE = EBP_boxcox_hcr_para_emprmse,
  ParaEstRMSE = EBP_boxcox_hcr_para_estrmse,
  WildEstRMSE = EBP_boxcox_hcr_wild_estrmse
) %>%
  pivot_longer(cols = -Domain, names_to = "Method", values_to = "RMSE")

p_hcr <- ggplot(df_hcr_long, aes(x = Domain, y = RMSE, color = Method)) +
  geom_line(size = 1) +
  geom_point() +
  labs(title = "HCR Indicator", x = "Domain", y = "RMSE") +
  scale_color_manual(values = c("EmpRMSE" = "black", 
                                "ParaEstRMSE" = "red", 
                                "WildEstRMSE" = "orange")) +
  theme_gray()

# ---- PGAP Indicator ----
df_pgap_long <- data.frame(
  Domain = Domain,
  EmpRMSE = EBP_boxcox_pgap_para_emprmse,
  ParaEstRMSE = EBP_boxcox_pgap_para_estrmse,
  WildEstRMSE = EBP_boxcox_pgap_wild_estrmse
) %>%
  pivot_longer(cols = -Domain, names_to = "Method", values_to = "RMSE")

p_pgap <- ggplot(df_pgap_long, aes(x = Domain, y = RMSE, color = Method)) +
  geom_line(size = 1) +
  geom_point() +
  labs(title = "PGAP Indicator", x = "Domain", y = "RMSE") +
  scale_color_manual(values = c("EmpRMSE" = "black", 
                                "ParaEstRMSE" = "red", 
                                "WildEstRMSE" = "orange")) +
  theme_gray()

# ---- QSR Indicator ----
df_qsr_long <- data.frame(
  Domain = Domain,
  EmpRMSE = EBP_boxcox_qsr_para_emprmse,
  ParaEstRMSE = EBP_boxcox_qsr_para_estrmse,
  WildEstRMSE = EBP_boxcox_qsr_wild_estrmse
) %>%
  pivot_longer(cols = -Domain, names_to = "Method", values_to = "RMSE")

p_qsr <- ggplot(df_qsr_long, aes(x = Domain, y = RMSE, color = Method)) +
  geom_line(size = 1) +
  geom_point() +
  labs(title = "QSR Indicator", x = "Domain", y = "RMSE") +
  scale_color_manual(values = c("EmpRMSE" = "black", 
                                "ParaEstRMSE" = "red", 
                                "WildEstRMSE" = "orange")) +
  theme_gray()+
  coord_cartesian(ylim = c(0, 4000))

# ---- Combine All Plots with Shared Legend ----
(p_mean + p_hcr + p_pgap + p_qsr) +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

