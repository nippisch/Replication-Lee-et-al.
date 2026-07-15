# JcAIC
library(MASS)
library(magic)
library(nlme)
library(emdi)
source("Lee_et_al_JcAIC/bias_correction_boot.R")
source("Lee_et_al_JcAIC/stepJcAIC.R")
source("Lee_et_al_JcAIC/optimal_parameter.R")

load("data/test_normal_sample.RData")

# Without Trafo-----------------------------------------------------------------


fixed <- as.formula(paste("y ~ ", paste(colnames(data)[3:8], collapse= "+")))

fitfull <- lme(fixed = fixed, random = ~ 1|domain, 
               method = "ML", data = data)
fitnull <- lme(fixed = y ~ 1, random = ~ 1|domain, 
               method = "ML", data = data)

B <- 5
interval <- c(-1, 2)
# object <- fitfull
# direction <- "both"
# scope <- list(lower = fitnull, upper = fitfull)
# trans <- "no"
# ycol <- "y"
# trace <- TRUE
# steps <- 1000

fitopt <- stepJcAIC(object = fitfull, direction = "both", 
                    scope = list(lower = fitnull, upper = fitfull),
                    B = B, trans = "no",
                    y_col = "y", 
                    interval = interval,
                    trace = TRUE)

summary(fitopt$fit)
JcAIC_no <- fitopt$models[[length(fitopt$models)]]$JcAIC

opt_fixed_no <- fitopt$fit$call$fixed
#-------------------------------------------------------------------------------

# Log Trafo---------------------------------------------------------------------
summary(data$y)
data$log_y <- log(data$y)
fixed_log <- update.formula(fixed, log_y ~.)

fitfull_log <- lme(fixed = fixed_log, random = ~ 1|domain, 
                   method = "ML", data = data)
fitnull_log <- lme(fixed = log_y ~ 1, random = ~ 1|domain, 
                   method = "ML", data = data)

fitopt_log <- stepJcAIC(object = fitfull_log, direction = "both", 
                        scope = list(lower = fitnull_log, upper = fitfull_log),
                        B = B, trans = "log",
                        y_col = "y", 
                        interval = interval,
                        trace = TRUE)

summary(fitopt_log$fit)
JcAIC_log <- fitopt_log$models[[length(fitopt_log$models)]]$JcAIC

opt_fixed_log <- fitopt_log$fit$call$fixed
#-------------------------------------------------------------------------------

# BoxCox Trafo------------------------------------------------------------------
summary(data$y)
shift <- 0
# optimal lambda for the full model
lambdaopt <- optimal_parameter(generic_opt,
                               fixed = fixed,
                               smp_data = data,
                               smp_domains = "domain",
                               transformation = "box.cox",
                               interval = interval)
# Initial box-cox transformed y with the full model
data$bc_y <- ((data$y + shift)^lambdaopt - 1)/lambdaopt
fixed_BC <- update.formula(fixed, bc_y ~.)

fitfull_BC <- lme(fixed = fixed_BC, random = ~ 1|domain, 
                  method = "ML", data = data)
fitnull_BC <- lme(fixed = bc_y ~ 1, random = ~ 1|domain, 
                  method = "ML", data = data)

fitopt_BC <- stepJcAIC(object = fitfull_BC, direction = "both", 
                       scope = list(lower = fitnull_BC, upper = fitfull_BC),
                       B = B, trans = "box.cox",
                       y_col = "y", 
                       interval = interval,
                       trace = TRUE)

summary(fitopt_BC$fit)
JcAIC_BC <- fitopt_BC$models[[length(fitopt_BC$models)]]$JcAIC

opt_fixed_BC <- fitopt_BC$fit$call$fixed

#-------------------------------------------------------------------------------
# Optimal model
c(opt_fixed_no, opt_fixed_log, opt_fixed_BC)

which.min(c(JcAIC_no, JcAIC_log, JcAIC_BC))