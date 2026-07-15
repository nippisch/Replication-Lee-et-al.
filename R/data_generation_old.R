###############################################################################
################### Data Generation Script ####################################
################### Author: Niklas Ippisch ####################################
###############################################################################

# loading relevant packages and functions --------------------------------------
library(dplyr)
source("Lee_et_al_JcAIC/bias_correction_boot.R")
source("Lee_et_al_JcAIC/stepJcAIC.R")
source("Lee_et_al_JcAIC/optimal_parameter.R")

# definition of global parameters ----------------------------------------------

N <- 10000 # population size
R <- 4 # monte Carlo iterations
D <- 50 # number of domains
B <- 5 # number of bootstrap iterations for JcAIC bias correction term
interval <- c(-2, 2) # interval for parameter lambda
base_seed <- 123456 # base seed for reproducibility 
dgp_settings <- c("normal_low", "normal_high", "log", "boxcox") # data generating processes
transformations <- c("no", "log", "boxcox_naive", "box.cox") # possible transformations
sim_res_sample <- rep(list(vector("list", R)), length(dgp_settings)) # results list for storing the samples

# results list for storing the JcAIC-results
sim_res_JcAIC <- setNames( 
  lapply(dgp_settings, function(d) {
    setNames(
      lapply(transformations, function(t) {
        vector("list", R)
      }),
      transformations
    )
  }),
  dgp_settings
)

# definition of help functions -------------------------------------------------

# sampling function
sampler <- function(DAT, sample_size, D) {
  
  # collecting all indexes of to sample observations
  idx <- unlist(lapply(1:D, function(d) {
    if (sample_size[d] == 0) return(NULL)  # skip out-of-sample domains
    rows <- which(DAT$domain == d)
    rows[sample(length(rows), size = min(sample_size[d], length(rows)))]
  }))
  
  # subsetting all to sample observations
  DAT[idx, ]
}

# data generation process function
dgp <- function(distribution, N, D) {
  
  # creating covariates independent of dgp
  domain <- sample(1:D, size = N, replace = TRUE)
  x3 <- rnorm(n = N, mean = 0, sd = 1)
  x4 <- rbeta(n = N, shape1 = 2, shape2 = 2)
  z1 <- rnorm(n = N, mean = runif(n = 1, min = 95, max = 105), sd = 20)
  z2 <- rbinom(n = N, size = 1, prob = 0.4)
  
  if (distribution == "normal_low") {
    x1 <- rnorm(n = N, mean = runif(n = 1, min = -3, max = 3), sd = 3)
    x2 <- rbinom(n = N, size = 1, prob = 0.2)
    
    u <- rnorm(n = D, mean = 0, sd = 30)
    e <- rnorm(n = N, mean = 0, sd = 60)
    
    y  <- 1200 - 10 * x1 + 100 * x2 - 10 * x3 + 25 * x4 + u[domain] + e
    
  } else if (distribution == "normal_high") {
    
    x1 <- rnorm(n = N, mean = runif(n = 1, min = -3, max = 3), sd = 3)
    x2 <- rbinom(n = N, size = 1, prob = 0.2)
    
    u <- rnorm(n = D, mean = 0, sd = 10)
    e <- rnorm(n = N, mean = 0, sd = 20)
    
    y  <- 1200 - 10 * x1 + 100 * x2 - 10 * x3 + 25 * x4 + u[domain] + e
    
  } else if (distribution == "log") {
    x1 <- rnorm(n = N, mean = runif(n = 1, min = 2, max = 3), sd = 2)
    x2 <- rbinom(n = N, size = 1, prob = 0.8)
    
    u <- rnorm(n = D, mean = 0, sd = 0.4)
    e <- rnorm(n = N, mean = 0, sd = 0.8)
    
    y  <- exp(10 - x1 + x2 - 0.5 * x3 - 0.3 * x4 + u[domain] + e)
    
  } else if (distribution == "boxcox") {
    
    x1 <- rnorm(n = N, mean = runif(n = 1, min = 2, max = 3), sd = 2)
    x2 <- rbinom(n = N, size = 1, prob = 0.8)
    
    u <- rnorm(n = D, mean = 0, sd = 0.4)
    e <- rnorm(n = N, mean = 0, sd = 0.8)
    
    y  <- ((10 - x1 + x2 - 0.5 * x3 - 0.3 * x4 + u[domain] + e) * (-0.5) + 1)^(1/-0.5)
    
  }
  
  data.frame(id = 1:N, domain, y, x1, x2, x3, x4, z1, z2)
}

# function for conducting the stepwise variable selection using JcAIC
stepwise_JcAIC <- function(data, B , interval, transformation = c("no", "log", "box.cox", "boxcox_naive")) {
  
  # building the base formula
  fixed <- as.formula(paste("y ~ ", paste(setdiff(colnames(data), c("id", "domain", "y")), collapse = "+")))
  
  if (transformation == "log") { # if statement 1: log-transformation
    
    # transformation and saving in new variable
    data$y_log <- log(data$y)
    
    # updating model formula
    fixed <- update.formula(fixed, y_log ~ .)
    
    # defining empty model
    fitnull_fixed <- as.formula("y_log ~ 1")
    
  } else if (transformation == "box.cox") { # if statement 2: boxcox-transformation
    
    # optimize lambda on full model using original y
    lambdaopt <- optimal_parameter(generic_opt,
                                   fixed = fixed,
                                   smp_data = data,
                                   smp_domains = "domain",
                                   transformation = "box.cox",
                                   interval = interval)
    
    # transformation and saving in new variable
    data$y_bc <- ((data$y)^lambdaopt - 1) / lambdaopt
    
    # updating model formula
    fixed <- update.formula(fixed, y_bc ~ .)
    
    # defining empty model
    fitnull_fixed <- as.formula("y_bc ~ 1")
    
  } else { # no transformation
    
    fitnull_fixed <- as.formula("y ~ 1")
    
  } # end if statement 2
  
  # define full model
  fitfull <- do.call(lme, list(fixed = fixed, 
                               random = ~ 1|domain, 
                               method = "ML", 
                               data = data))
  
  # define empty model
  fitnull <- do.call(lme, list(fixed = fitnull_fixed, 
                               random = ~ 1|domain, 
                               method = "ML", 
                               data = data))
  
  trans_tmp <- trans_tmp <- if (transformation == "boxcox_naive") "no" else transformation
  
  # conduct stepwise variable selection
  fitopt <- stepJcAIC(object = fitfull, 
                      direction = "both", 
                      scope = list(lower = fitnull, upper = fitfull),
                      B = B, 
                      trans = trans_tmp,
                      y_col = "y",
                      interval = interval,
                      trace = TRUE)
  
  if (transformation != "boxcox_naive") {
    
    JcAIC_res <- fitopt$models[[length(fitopt$models)]]
    JcAIC_res$opt_vars <- fitopt$fit$call$fixed
    
  } else if (transformation == "boxcox_naive") {
    
    # optimize lambda on full model using original y
    lambdaopt <- optimal_parameter(generic_opt,
                                   fixed = as.formula(fitopt$fit$call$fixed),
                                   smp_data = data,
                                   smp_domains = "domain",
                                   transformation = "box.cox",
                                   interval = interval)
    
    # transformation and saving in new variable
    data$y_bc_naive <- ((data$y)^lambdaopt - 1) / lambdaopt
    
    fixed_updated <- update.formula(as.formula(fitopt$fit$call$fixed), y_bc_naive ~ .)
    
    fit_updated <- do.call(lme, list(fixed = fixed_updated, 
                                 random = ~ 1|domain, 
                                 method = "ML", 
                                 data = data))
    
    JcAIC_result <- JcAIC(lmefit = fit_updated, domain = "domain", trans = "box.cox", m = NA, 
                          interval = interval, B = B, y_col = "y")
    
    JcAIC_res <- list(lambdaopt = lambdaopt,
                      JcAIC = JcAIC_result$JcAIC,
                      JcAIC_result = JcAIC_result,
                      opt_vars = fixed_updated)
    
  }
  
  # return results
  return(JcAIC_res)
  
}

# simulation study -------------------------------------------------------------

for (r in 1:R) { # loop level 1: R monte carlo iterations
  
  for (j in 1:length(dgp_settings)) { # loop level 2: data generating processes
    
    # varying but reproducible seed
    set.seed(base_seed + 1000 * r + j)
    
    # generating the population
    pop_r <- dgp(distribution = dgp_settings[j], N = N, D = D)
    
    # defining the sample size, drawing the sample, and saving it
    raw_sizes <- round(runif(D, min = 0, max = 30))
    sample_size <- ifelse(raw_sizes < 3, 0, pmax(raw_sizes, 5)) # prevent areas with just 1-2 observations to avoid non-convergence of LMM
    sample_r <- sampler(DAT = pop_r, sample_size = sample_size, D = D)
    sim_res_sample[[j]][[r]] <- sample_r
    
    # fitting LMM and variable selection
    
    for (t in 1:length(transformations)) { # loop level 3: iterating over transformations
      
      res <- stepwise_JcAIC(data = sample_r, B = B, interval = interval, transformation = transformations[t])
      
      sim_res_JcAIC[[dgp_settings[j]]][[transformations[t]]][[r]] <- res
      
      cat(paste("Finished: R = ", r, " DGP = ", dgp_settings[j], " Transformation: ", transformations[t], "\n"))
      
    } # end loop level 3
    
    
  } # end loop level 2
  
} # end loop level 1

save(sim_res_JcAIC, sim_res_sample, sample_r, pop_r, file = "data/test_r4_final.RData")

# Vergleich mit direktem Schätzer in areas mit großer Stichprobe
# Vergleich von small area estimates in aggreggierter Form --> auf höherer Ebene guter design-basierter Schätzer
# Unsicherheit: CV's von Schätzern
# Dezidiert in Hausarbeit, wie funktioniert Vergleich ohne bekannte Wahrheit

