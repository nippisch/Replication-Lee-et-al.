# loading relevant packages and functions --------------------------------------
library(future)
library(future.apply)
library(progressr)
library(nlme)
library(dplyr)
source("Lee_et_al_JcAIC/bias_correction_boot.R")
source("Lee_et_al_JcAIC/stepJcAIC.R")
source("Lee_et_al_JcAIC/optimal_parameter.R")

# definition of global parameters ----------------------------------------------

N <- 10000 # population size
R <- 6 # monte Carlo iterations
D <- 50 # number of domains
B <- 50 # number of bootstrap iterations for JcAIC bias correction term
interval <- c(-2, 2) # interval for parameter lambda
base_seed <- 12345 # base seed for reproducibility 
set.seed(base_seed) # set base seed for reproducibility
sample_size <- round(runif(D, min = 0, max = 30)) # vector of sample sizes per domain
dgp_settings <- c("normal_low", "normal_high") # data generating processes
transformations <- c("box.cox") # possible transformations

# definition of help functions -------------------------------------------------

# sampling function
sampler <- function(data, sample_size, D) {
  
  # collecting all indexes of to sample observations
  idx <- unlist(lapply(1:D, function(d) {
    rows <- which(data$domain == d)
    rows[sample(length(rows), size = min(sample_size[d], length(rows)))]
  }))
  
  # subsetting all to sample observations
  data[idx, ]
}

# data generation process function
dgp <- function(distribution, N, D) {
  
  # creating covariates independent of dgp
  domain <- sample(1:D, size = N, replace = TRUE)
  x2 <- rbinom(n = N, size = 1, prob = 0.8)
  x3 <- rnorm(n = N, mean = 0, sd = 1)
  z1 <- rnorm(n = N, mean = 1, sd = 0.1)
  z2 <- rbinom(n = N, size = 1, prob = 0.4)
  
  if (distribution == "normal_low") {
    
    x1 <- rnorm(n = N, mean = runif(n = 1, min = -3, max = 3), sd = 3)
    
    u <- rnorm(n = D, mean = 0, sd = 30)
    e <- rnorm(n = N, mean = 0, sd = 60)
    
    y  <- 400 - 10 * x1 + 100 * x2 - 10 * x3 + u[domain] + e
    
  } else if (distribution == "normal_high") {
    
    x1 <- rnorm(n = N, mean = runif(n = 1, min = -3, max = 3), sd = 3)
    
    u <- rnorm(n = D, mean = 0, sd = 10)
    e <- rnorm(n = N, mean = 0, sd = 20)
    
    y  <- 400 - 10 * x1 + 100 * x2 - 10 * x3 + u[domain] + e
    
  } else if (distribution == "log") {
    
    x1 <- rnorm(n = N, mean = runif(n = 1, min = 2, max = 3), sd = 2)
    
    u <- rnorm(n = D, mean = 0, sd = 0.4)
    e <- rnorm(n = N, mean = 0, sd = 0.8)
    
    y  <- exp(10 - x1 + x2 - 0.5 * x3 + u[domain] + e)
    
  } else if (distribution == "boxcox") {
    
    x1 <- rnorm(n = N, mean = runif(n = 1, min = 2, max = 3), sd = 2)
    
    u <- rnorm(n = D, mean = 0, sd = 0.4)
    e <- rnorm(n = N, mean = 0, sd = 0.8)
    
    y  <- ((10 - x1 + x2 - 0.5 * x3 + u[domain] + e) * (-0.5) + 1)^(1/-0.5)
    
  }
  
  data.frame(id = 1:N, domain, y, x1, x2, x3, z1, z2)
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
    lambdaopt_start <- optimal_parameter(generic_opt,
                                         fixed = fixed,
                                         smp_data = data,
                                         smp_domains = "domain",
                                         transformation = "box.cox",
                                         interval = interval)
    
    # transformation and saving in new variable
    data$y_bc <- ((data$y)^lambdaopt_start - 1) / lambdaopt_start
    
    # updating model formula
    fixed <- update.formula(fixed, y_bc ~ .)
    
    # defining empty model
    fitnull_fixed <- as.formula("y_bc ~ 1")
    
  } else { # no transformation
    
    fitnull_fixed <- as.formula("y ~ 1")
    
  } # end if statement 2
  
  # loosen optimizer controls to avoid non-convergence
  ctrl <- lmeControl(maxIter = 200, msMaxIter = 200,
                     niterEM = 100, msMaxEval = 400, apVar = FALSE)
  
  # define full model
  fitfull <- do.call(lme, list(fixed = fixed, 
                               random = ~ 1|domain, 
                               method = "ML", 
                               data = data,
                               control = ctrl))
  
  # define empty model
  fitnull <- do.call(lme, list(fixed = fitnull_fixed, 
                               random = ~ 1|domain, 
                               method = "ML", 
                               data = data,
                               control = ctrl))
  
  # setting transformation on "no" for boxcox_naive for the stepJcAIC-Function
  trans_tmp <- if (transformation == "boxcox_naive") "no" else transformation
  
  # conduct stepwise variable selection
  fitopt <- stepJcAIC(object = fitfull, 
                      direction = "both", 
                      scope = list(lower = fitnull, upper = fitfull),
                      B = B, 
                      trans = trans_tmp,
                      y_col = "y",
                      interval = interval,
                      trace = FALSE)
  
  if (transformation != "boxcox_naive") {
    
    JcAIC_res <- fitopt$models[[length(fitopt$models)]]
    JcAIC_res$opt_vars <- fitopt$fit$call$fixed
    JcAIC_res$lambdaopt_start <- lambdaopt_start
    
  } else if (transformation == "boxcox_naive") {
    
    # optimize lambda on optimal model
    lambdaopt_start <- optimal_parameter(generic_opt,
                                         fixed = formula(fitopt$fit),
                                         smp_data = data,
                                         smp_domains = "domain",
                                         transformation = "box.cox",
                                         interval = interval)
    
    # transformation and saving in new variable
    data$y_bc_naive <- ((data$y)^lambdaopt_start - 1) / lambdaopt_start
    
    fixed_updated <- update.formula(formula(fitopt$fit), y_bc_naive ~ .)
    
    fit_updated <- do.call(lme, list(fixed = fixed_updated, 
                                     random = ~ 1|domain, 
                                     method = "ML", 
                                     data = data))
    
    # calculating final JcAIC for boxcox_naive setup
    JcAIC_result <- JcAIC(lmefit = fit_updated, 
                          domain = "domain", 
                          trans = "box.cox", 
                          m = NA, 
                          interval = interval, 
                          B = B, 
                          y_col = "y")
    
    JcAIC_res <- list(lambdaopt = lambdaopt_start,
                      lambdaopt_start = lambdaopt_start,
                      JcAIC = JcAIC_result$JcAIC,
                      JcAIC_result = JcAIC_result,
                      opt_vars = fixed_updated)
    
  }
  
  # return results
  return(JcAIC_res)
  
}

# function for one iteration
run_one_iteration <- function(r, dgp_settings, transformations, N, D, B,
                              interval, base_seed, sample_size, p = NULL) {
  
  # Pre-allocate output for this iteration
  iter_sample <- setNames(vector("list", length(dgp_settings)), dgp_settings)
  iter_JcAIC  <- setNames(
    lapply(dgp_settings, function(d) setNames(vector("list", length(transformations)), 
                                              transformations)),
    dgp_settings)
  
  for (j in seq_along(dgp_settings)) {
    
    # setting seed for reproducibility
    set.seed(base_seed + 1000 * r + j)
    
    # generating population and draw sample
    pop_r    <- dgp(distribution = dgp_settings[j], N = N, D = D)
    sample_r <- sampler(data = pop_r, sample_size = sample_size, D = D)
    iter_sample[[dgp_settings[j]]] <- sample_r
    
    for (t in seq_along(transformations)) {
      iter_JcAIC[[dgp_settings[j]]][[transformations[t]]] <- tryCatch(
        stepwise_JcAIC(data = sample_r, B = B, interval = interval, transformation = transformations[t]),
        error = function(e) {
          list(error = conditionMessage(e),
               r = r, dgp = dgp_settings[j], transformation = transformations[t])
        }
      )
    }
    
    # tracking process
    if (!is.null(p)) p(message = paste0("R = ", r, ", DGP = ", dgp_settings[j]))
  }
  
  list(sample = iter_sample, JcAIC = iter_JcAIC)
}

# parallel execution -----------------------------------------------------------

# initializing the setup
plan(multisession, workers = availableCores() - 2)
handlers(global = TRUE)
handlers("txtprogressbar")

# execute parallel
all_results <- with_progress({
  
  p <- progressor(steps = R * length(dgp_settings))
  
  future_lapply(
    X               = seq_len(R),
    FUN             = run_one_iteration,
    # global parameters
    dgp_settings    = dgp_settings,
    transformations = transformations,
    N               = N,
    D               = D,
    B               = B,
    interval        = interval,
    base_seed       = base_seed,
    sample_size     = sample_size,
    p               = p,
    # seed handling
    future.seed     = FALSE,
    # help functions
    future.globals  = c(
      "dgp", "sampler", "stepwise_JcAIC", "log_shift_std",
      "stepJcAIC", "JcAIC", "optimal_parameter", "generic_opt", "reml",
      "box_cox", "box_cox_std", "box_cox_back", "std_data_transformation",
      "boot_biascorrection_DD", "boot_biascorrection_DD_normal", 
      "boot_biascorrection_DD_log", "geometric.mean", "addterm_JcAIC",
      "bound_estimation", "conditional_loglikelihood_boot", "dropterm_JcAIC",
      "dual_power", "dual_power_std", "log_Jacobian_boxcox", "log_shift"),
    # necessary packages
    future.packages = c("nlme", "dplyr", "MASS", "magic", "Ecfun")
  )
  
})

# Reset to sequential execution when done
plan(sequential)

# store results in proper structure
sim_res_sample <- setNames(
  lapply(dgp_settings, function(d) lapply(seq_len(R), function(r) all_results[[r]]$sample[[d]])),
  dgp_settings
)

sim_res_JcAIC <- setNames(
  lapply(dgp_settings, function(d) {
    setNames(
      lapply(transformations, function(t) {
        lapply(seq_len(R), function(r) all_results[[r]]$JcAIC[[d]][[t]])
      }),
      transformations
    )
  }),
  dgp_settings
)

# storing parameters used
names(sample_size) <- 1:D
params <- list(N = N, D = D, R = R, interval_lambda = interval, sample_sizes = sample_size, dgps = dgp_settings, transformation = transformations, B = B)

# save results
save(params, sim_res_sample, sim_res_JcAIC, file = "data/simulation_study_results_raw.RData")
