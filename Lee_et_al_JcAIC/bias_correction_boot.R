library(Ecfun)
boot_biascorrection_DD_normal <- function(lmefit, Z, X, D, B){
  
  N <- nobs(lmefit)
  # step 1: get thetahat
  # variance of random effect, u
  sigma_u2_hat <- as.numeric(nlme::VarCorr(lmefit)[1,1])
  # variance of individual error, e
  sigma_e2_hat <- lmefit$sigma^2
  
  B1 <- B2 <- c()
#start_time <- Sys.time()  
  for (b in 1:B){
    set.seed(b)
    # b = 1
    # setp 2: generate random effect, u and individual error e
    e_boot <- rnorm(N , 0, sqrt(sigma_e2_hat))
    u_boot <- rnorm(D , 0, sqrt(sigma_u2_hat))
    
    # step 3: bootstrap model 
    
    y_boot <- X %*% fixed.effects(lmefit) + Z %*% u_boot + e_boot
    

    # step 4: estimate bootstrap estimators: theta_boot
    dat_boot <- lmefit$data

    #fixed <- as.formula(as.character(lmefit$call$fixed)) # Yeonjoo: replaced by 28
    fixed <- formula(lmefit$call$fixed)
    dat_boot[[fixed[[2]]]] <- y_boot
    
    environment(fixed) <- environment()
    lmefit_boot <- lme(fixed = fixed, random = as.formula(lmefit$call$random), 
                       method = "ML", data = dat_boot, 
                       control = lmeControl(opt = "optim"))
    lmefit_boot$call$fixed <- fixed
    lmefit_boot$call$random <- as.formula(lmefit$call$random)
 ###############################################################################   
    clogLik_boot_diff <- conditional_loglikelihood_boot(lmefit = lmefit_boot,
                                              y_boot = lmefit$data[[fixed[[2]]]],
                                              y_starboot = as.vector(dat_boot[[fixed[[2]]]]),
                                              dat_boot = lmefit$data,
                                              dat_starboot = dat_boot)
################################################################################
    # condLogLik of y_tilde given theta_boot
    #clogLik_boot <- conditional_loglikelihood(lmefit = lmefit_boot,
    #                                          y = lmefit$data[[fixed[[2]]]],
    #                                          dat = lmefit$data)
    # condLogLik of y_tilde_boot given theta_boot
    #clogLik_starboot <- conditional_loglikelihood(lmefit = lmefit_boot,
    #                                              y = as.vector(dat_boot[[fixed[[2]]]]),
    #                                              dat = dat_boot)
    
    
    B1[b] <- clogLik_boot_diff
    #B1[b] <- (clogLik_starboot) - (clogLik_boot) 
  }
#end_time <- Sys.time()
#time_dur <- end_time - start_time
#time_dur

  list(B1 = mean(B1),
       B2 = NA)
}

boot_biascorrection_DD_log <- function(lmefit, Z, X, D, B, m, y_col){
  
  N <- nobs(lmefit)
  # step 1: get thetahat
  # variance of random effect, u
  sigma_u2_hat <- as.numeric(nlme::VarCorr(lmefit)[1,1])
  # variance of individual error, e
  sigma_e2_hat <- lmefit$sigma^2
  
  B1 <- B2 <- c()
  for (b in 1:B){
    set.seed(b)
    # b = 1
    # setp 2: generate random effect, u and individual error e
    e_boot <- rnorm(N , 0, sqrt(sigma_e2_hat))
    u_boot <- rnorm(D , 0, sqrt(sigma_u2_hat))
    
    # step 3: bootstrap model 
    
    logy_boot <- X %*% fixed.effects(lmefit) + Z %*% u_boot + e_boot
    
    
    # step 4: estimate bootstrap estimators: theta_boot
    dat_boot <- lmefit$data
    
    fixed <- formula(lmefit$call$fixed)
    dat_boot[[fixed[[2]]]] <- as.numeric(logy_boot)
    dat_boot$y <- exp(dat_boot[[fixed[[2]]]])
    
    environment(fixed) <- environment()
    lmefit_boot <- lme(fixed = fixed, random = as.formula(lmefit$call$random), 
                       method = "ML", data = dat_boot, control = lmeControl(opt = "optim"))
    lmefit_boot$call$fixed <- fixed
    lmefit_boot$call$random <- as.formula(lmefit$call$random)
    
    ###############################################################################   
    clogLik_boot_diff <- conditional_loglikelihood_boot(lmefit = lmefit_boot,
                                                        y_boot = lmefit$data[[fixed[[2]]]],
                                                        y_starboot = as.vector(dat_boot[[fixed[[2]]]]),
                                                        dat_boot = lmefit$data,
                                                        dat_starboot = dat_boot)
    ################################################################################
    
    # condLogLik of y_tilde given theta_boot
    #clogLik_boot <- conditional_loglikelihood(lmefit = lmefit_boot,
    #                                          y = lmefit$data[[fixed[[2]]]],
    #                                          dat = lmefit$data)
    # condLogLik of y_tilde_boot given theta_boot
    #clogLik_starboot <- conditional_loglikelihood(lmefit = lmefit_boot,
    #                                              y = as.vector(dat_boot[[fixed[[2]]]]),
    #                                              dat = dat_boot)
    
    
    B1[b] <- clogLik_boot_diff - sum(log(dat_boot$y + m)) + 
      sum(log(lmefit$data[[y_col]] + m)) 

    #B1[b] <- (clogLik_starboot - sum(log(dat_boot$y + m))) - 
      #(clogLik_boot - sum(log(lmefit$data$y + m))) 
   
  }
  
  list(B1 = mean(B1),
       B2 = NA)
}

boot_biascorrection_DD <- function(lmefit, Z, X, D, lambdaopt, m, interval = c(-2, 2), B, y_col){
  
  N <- nobs(lmefit)
  # step 1: get thetahat
  # variance of random effect, u
  sigma_u2_hat <- as.numeric(nlme::VarCorr(lmefit)[1,1])
  # variance of individual error, e
  sigma_e2_hat <- lmefit$sigma^2
  
  B1 <- B2 <- c()
  lambdaB <- c()
  for (b in 1:B){
    set.seed(b)
    #print(b)
    # b = 1
    # setp 2: generate random effect, u and individual error e
    e_boot <- rnorm(N , 0, sqrt(sigma_e2_hat))
    u_boot <- rnorm(D , 0, sqrt(sigma_u2_hat))
    
    # step 3: bootstrap model and backtransform
    
    y_bc_boot <- X %*% nlme::fixed.effects(lmefit) + Z %*% u_boot + e_boot
    
    #y_boot <- box_cox_back(y = y_bc_boot, lambda = lambdaopt)[, 1]
    # gibt NA: Problem wie hier
    # https://stats.stackexchange.com/questions/541748/simple-problem-with-box-cox-transformation-in-a-time-series-model
    
    #y_boot <- Ecfun::invBoxCox(y_bc_boot, lambdaopt)
    y_boot <- (as.complex(y_bc_boot) * lambdaopt + 1)^(1/lambdaopt) - m
    
    # step 4: estimate bootstrap estimators: theta_boot
    dat_boot <- lmefit$data
    
    dat_boot[[lmefit$call$fixed[[2]]]] <- y_bc_boot
    fixed <- formula(lmefit$call$fixed)
    
    lmefit_boot <- NULL
    try(lmefit_boot <- nlme::lme(fixed = fixed, random = as.formula(lmefit$call$random), 
                                 method = "ML", data = dat_boot), silent = TRUE)
    if (is.null(lmefit_boot)){
      lmefit_boot <- nlme::lme(fixed = fixed, 
                               random = as.formula(lmefit$call$random), 
                               method = "ML", data = dat_boot, 
                               control = nlme::lmeControl(opt = "optim"))
    } 
    lmefit_boot <- lmefit_boot
    
    #dat_boot$y <- y_boot
    dat_starboot <- lmefit$data
    dat_starboot$y <- as.numeric(y_boot)
    
    #fixed <- as.formula(as.character(lmefit$call$fixed)) 
    fixed <- formula(lmefit$call$fixed) 
    environment(fixed) <- environment()
    
    #lamdaopt_starboot <- NULL
    #try(
    lamdaopt_starboot <- optimal_parameter(generic_opt,
                                           fixed = as.formula(update.formula(fixed,
                                                                             as.formula(paste("y", "~ .")))),
                                           smp_data = dat_starboot,
                                           smp_domains = names(lmefit$coefficients$random),
                                           transformation = "box.cox",
                                           interval = interval)
    #, silent = TRUE)
    
    
    fixed <- formula(lmefit$call$fixed)
    environment(fixed) <- environment()
    
    dat_starboot[[fixed[[2]]]] <-  box_cox(lamdaopt_starboot, dat_starboot$y)$y
    m_starboot <- box_cox(lamdaopt_starboot, dat_starboot$y)$m
    
    lmefit_starboot <- NULL
    try(lmefit_starboot <- nlme::lme(fixed = fixed, random = as.formula(lmefit$call$random), 
                                     method = "ML", data = dat_starboot), silent = TRUE)
    if (is.null(lmefit_starboot)){
      lmefit_starboot <- nlme::lme(fixed = fixed, 
                                   random = as.formula(lmefit$call$random), 
                                   method = "ML", data = dat_starboot, 
                                   control = nlme::lmeControl(opt = "optim"))
    } 
    lmefit_starboot <- lmefit_starboot
    
    lmefit_boot$call$fixed <- formula(lmefit$call$fixed)
    lmefit_boot$call$random <- as.formula(lmefit$call$random)
    lmefit_starboot$call$fixed <- formula(lmefit$call$fixed)
    lmefit_starboot$call$random <- as.formula(lmefit$call$random)
    
    ###############################################################################   
    clogLik_boot_diff <- conditional_loglikelihood_boot(lmefit = lmefit_boot,
                                                        y_boot = lmefit$data[[fixed[[2]]]],
                                                        y_starboot = as.vector(dat_boot[[fixed[[2]]]]),
                                                        dat_boot = lmefit$data,
                                                        dat_starboot = dat_boot)
    ################################################################################
    
    # condLogLik of y_tilde given theta_boot
    #clogLik_boot <- conditional_loglikelihood(lmefit = lmefit_boot,
    #                                          y = lmefit$data[[fixed[[2]]]],
    #                                          dat = lmefit$data)
    
    # condLogLik of y_tilde_boot given theta_boot
    #clogLik_starboot <- conditional_loglikelihood(lmefit = lmefit_starboot,
    #                                              y = as.vector(dat_starboot[[fixed[[2]]]]),
    #                                              dat = dat_starboot)
    
    
    logJacobian_boot <- log_Jacobian_boxcox(lambdaopt, y = lmefit$data[[y_col]], m = m)
    logJacobian_starboot <- log_Jacobian_boxcox(lamdaopt_starboot, y = dat_starboot$y, m = m_starboot)
    
    B1[b] <- clogLik_boot_diff + logJacobian_starboot - logJacobian_boot
    #B1[b] <- (clogLik_starboot + logJacobian_starboot) - (clogLik_boot + logJacobian_boot)
    #B1[b] <- (clogLik_boot + logJacobian_boot) - (clogLik_starboot + logJacobian_starboot) # Yeonjoo - test
    lambdaB[b] <- lamdaopt_starboot
  }
  
  list(B1 = mean(B1),
       B2 = NA)
}

conditional_loglikelihood_boot <- function(lmefit,
                                           y_boot,
                                           y_starboot,
                                           dat_boot,
                                           dat_starboot){
  # Conditional loglikelihood
  D <- length(unique(lmefit$data[[names(lmefit$coefficients$random)]]))
  N <- nobs(lmefit)
  
  # variance of random effect, u
  sigma_u2 <- as.numeric(nlme::VarCorr(lmefit)[1,1])
  # variance of individual error, e
  sigma_e2 <- lmefit$sigma^2
  
  inv_conV_y2 <- diag((sigma_e2)^(-1), ncol = N, nrow = N)
  
  mu_boot <- as.numeric(nlme:::predict.lme(lmefit, newdata = dat_boot))
  mu_starboot <- as.numeric(nlme:::predict.lme(lmefit, newdata = dat_starboot))
  
  clogLik_boot <- -0.5*(N*log(2*pi) + N*log(sigma_e2) + t(y_boot - mu_boot) %*%
                          inv_conV_y2 %*% (y_boot - mu_boot))
  
  clogLik_starboot <- -0.5*(N*log(2*pi) + N*log(sigma_e2) + t(y_starboot - mu_starboot) %*%
                          inv_conV_y2 %*% (y_starboot - mu_starboot))
  
  diff <- (clogLik_starboot) - (clogLik_boot)
  diff
}

log_Jacobian_boxcox <- function(lambdaopt, y, m){
  (lambdaopt - 1) * sum(log(y + m))
}