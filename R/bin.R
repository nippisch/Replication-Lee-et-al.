###############################################################################
################### Data Generation Script ####################################
################### Author: Niklas Ippisch ####################################
###############################################################################

# loading relevant packages ----------------------------------------------------
library(tidyverse)
library(lme4)
library(emdi)
library(saeSim)

# definition of global parameters ----------------------------------------------


set.seed(123456)
R <- 10 # Monte Carlo iterations
D <- 50  # Number of Small Areas
sample_size <- rpois(R, lambda = 22) # area level specific sample sizes
summary(sample_size)
N <- rep(200,Domains)

dgp <- c("normal", "log", "pareto") # data generating processes


# helper functions for covariates ----------------------------------------------

gen_x1 <- function(dat) {
  dat["x1"] <- rnorm(n = N, mean = runif(n = 1, min = -3, max = 3), sd = 3^2)
  dat
}

gen_x2 <- function(dat) {
  dat["x2"] <- rbinom(n = N, size = 1, prob = 0.8)
  dat
}

gen_x3 <- function(dat) {
  dat["x3"] <- rpois(n = N, lambda = 15)
  dat
}

gen_x4 <- function(dat) {
  dat["x4"] <- rnorm(n = N, mean = 0, sd = 5^2)
  dat
}

gen_z1 <- function(dat) {
  dat["z1"] <- rnorm(n = N, mean = runif(n = 1, min = 95, max = 105), sd = 2^2)
  dat
}

gen_z2 <- function(dat) {
  dat["z2"] <- rbinom(n = N, size = 1, prob = 0.4)
  dat
}

# Define the components of one arbitrary population
setup <- sim_base(data = base_id(nDomains = D, nUnits = N)) |> 
  sim_gen(gen_generic(rnorm, mean = runif(n = 1, min = -3, max = 3), sd = 3^2, name = "x1", groupVars = "idD"))
sim_gen(gen_generic(runif, min = -3, max = 3, groupVars="idD", name = "muD")) |> 
  sim_gen(gen_XNorm) |> 
  as.data.frame |> 
  sim_gen(generator=gen_myE) |> 
  sim_gen_v(mean=0,sd = 500) |> 
  sim_resp_eq(y = 12000 - 400 * x + v + e)  

# Generate "simruns" populations defined by the setup
pop <- sim(setup, R = R)
# Set the negative values to 0 (to avoid negative income data)
Pop <- lapply(Pop, function(X){ X$y[X$y<0] = 0; return (X) })
# Generate "simruns" samples from the population with "sample_size" area-specific sample sizes
Pop <- lapply(Pop,sampler)




# Help functions for x, e and sampling
gen_XNorm=function(dat,m=dat$muD, s=7.5){
  dat["x"]=rnorm(nrow(dat), mean = m, sd = s)
  return(dat)
}
gen_myE=function(dat,shape=3, scale=2000){
  tmp =  sqrt(2) * rpareto(nrow(dat), shape = shape, scale = scale)
  dat["e"]= tmp - mean(tmp)
  return(dat)
}
gen_myE <- function(dat,m=0, s=1000) {
  dat["e"]=rnorm(nrow(dat), mean = m, sd = s)
  dat
}
sampler <- function(DAT){
  smp=as.data.frame(matrix(nrow=sum(sample_size) , ncol=ncol(DAT)))
  brd=append(0,cumsum(sample_size))
  for(i in 1:Domains){
    smp[((brd[i]+1):brd[i+1]),]=(DAT[DAT$idD==i,])[sample(1:sum(DAT$idD==i),size=sample_size[i]),]
  }
  attr(smp,"pop")=DAT
  colnames(smp)=colnames(DAT)
  return(smp)
}



# simulation study -------------------------------------------------------------

for (i in 1:R) { # loop level 1: R monte carlo iterations
  
  set.seed(base_seed + i)
  
  # creating the covariates
  x1 <- rnorm(n = N, mean = runif(n = 1, min = -3, max = 3), sd = 9)
  x2 <- rbinom(n = N, size = 1, prob = 0.8)
  x3 <- rpois(n = N, lambda = 15)
  x4 <- rnorm(n = N, mean = 0, sd = 5^2)
  z1 <- rnorm(n = N, mean = runif(n = 1, min = 95, max = 105), sd = 5)
  z2 <- rbinom(n = N, size = 1, prob = 0.4)
  
  for (j in 1:length(dgp)) { # loop level 2: 3 data generating processes
    
    if (dgp[j] == "normal") { # if statement 1: normal distribution
      
      # randomly putting observations in one of the D domains
      domain <- sample(1:D, size = N, replace = TRUE)
      
      u <- rnorm(D, mean = 0, sd = 10)
      e <- rnorm(N, mean = 0, sd = 20)
      
      # creating dependent variable and binding it into data frame
      y <- 500 - 10*x1 + 100*x2 + 30*x3 + 0.03*x4 + u[domain] + e
      
      df_r <- data.frame(y, x1, x2, x3, x4, z1, z2, domain = domain)
      
      # creating sample sizes proportional to domain size
      domain_sizes <- table(df_r$domain)
      domain_sizes <- as.numeric(domain_sizes[as.character(1:D)])
      p_h <- pmin(domain_sizes / sum(domain_sizes) * 550 / 30, 1) # 550 being approximately the aimed total sample size
      sample_sizes <- rbinom(D, size = 30, prob = p_h)
      sample_sizes <- pmin(sample_sizes, domain_sizes)
      
      # drawing the sample per domain
      df_r$sampled <- 0
      for (d in 1:D) { # loop level 3: iterating over rows with domain d
        domain_rows <- which(df_r$domain == d)
        if (sample_sizes[d] > 0) {
          chosen <- sample(domain_rows, size = sample_sizes[d])
          df_r$sampled[chosen] <- 1
        }
      }
      
      # Final sample
      sample_r <- df_r[df_r$sampled == 1, ]
      
      lmm_r <- lmer(y ~ x1 + x2 + x3 + x4 + z1 + z2 + (1 | domain), data = sample_r)
      
      
    } # end if statement 1
    
  } # end loop level 2
  
  
} # end loop level 1



# creating sample sizes proportional to domain size
domain_sizes <- table(df_r$domain)
domain_sizes <- as.numeric(domain_sizes[as.character(1:D)])
p_h <- pmin(domain_sizes / sum(domain_sizes) * 550 / 30, 1) # 550 being approximately the aimed total sample size
sample_sizes <- rbinom(D, size = 30, prob = p_h)
sample_sizes <- pmin(sample_sizes, domain_sizes)

# drawing the sample per domain
df_r$sampled <- 0
for (d in 1:D) { # loop level 3: iterating over rows with domain d
  domain_rows <- which(df_r$domain == d)
  if (sample_sizes[d] > 0) {
    chosen <- sample(domain_rows, size = sample_sizes[d])
    df_r$sampled[chosen] <- 1
  }
}

# Final sample
sample_r <- df_r[df_r$sampled == 1, ]

lmm_r <- lmer(y ~ x1 + x2 + x3 + x4 + z1 + z2 + (1 | domain), data = sample_r)

