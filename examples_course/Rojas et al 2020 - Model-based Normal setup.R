#########################################################################
# Advanced Small Area Estimation 
#
# Author:    Prof. Dr. Timo Schmid
# Content:   Data generating process for setting normal (Rojas-Perilla et al. 2020)
#########################################################################
rm(list=ls())

### Load required packages
library(emdi)
library(saeSim)

### Step 1: Data generation with R package -------------------------------------------------

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
  attr(dat, "pop")=dat
  dat
}

### b) Define the data generation mechanism -------------------------------------------------
# Normal setup compare Rojas-Perilla et al. (2020) Table 6 
# Help functions for x,e and sampling
gen_XNorm <- function(dat,m=dat$muD, s=3){
  dat["x"]=rnorm(nrow(dat), mean = m, sd = s)
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


# save(Pop, file = "Setting1_Normal_500.RData")