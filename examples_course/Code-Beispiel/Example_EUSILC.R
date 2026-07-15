library(MASS)
library(magic)
library(nlme)
library(emdi)

source("JcAIC_code/bias_correction_boot.R")
source("JcAIC_code/stepJcAIC.R")
source("JcAIC_code/optimal_parameter.R")

################################################################################
#                                   Data                                       #
################################################################################
data("eusilcA_smp")
summary(eusilcA_smp)
################################################################################
#                           optimal Model by JcAIC                              #
################################################################################
B <- 5
interval <- c(-1, 2)

# Without Trafo-----------------------------------------------------------------
colnames(eusilcA_smp)
fixed <- as.formula(paste("eqIncome ~ ", paste(colnames(eusilcA_smp)[2:14], collapse= "+")))

fitfull <- lme(fixed = fixed, random = ~ 1|district, 
               method = "ML", data = eusilcA_smp)
fitnull <- lme(fixed = eqIncome ~ 1, random = ~ 1|district, 
               method = "ML", data = eusilcA_smp)

fitopt <- stepJcAIC(object = fitfull, direction = "both", 
                    scope = list(lower = fitnull, upper = fitfull),
                    B = B, trans = "no",
                    y_col = "eqIncome", 
                    interval = interval,
                    trace = TRUE)

summary(fitopt$fit)
JcAIC_no <- fitopt$models[[length(fitopt$models)]]$JcAIC

opt_fixed_no <- fitopt$fit$call$fixed
#-------------------------------------------------------------------------------

# Log Trafo---------------------------------------------------------------------
summary(eusilcA_smp$eqIncome)
eusilcA_smp$log_eqIncome <- log(eusilcA_smp$eqIncome)
fixed_log <- update.formula(fixed, log_eqIncome ~.)

fitfull_log <- lme(fixed = fixed_log, random = ~ 1|district, 
                   method = "ML", data = eusilcA_smp)
fitnull_log <- lme(fixed = log_eqIncome ~ 1, random = ~ 1|district, 
                   method = "ML", data = eusilcA_smp)

fitopt_log <- stepJcAIC(object = fitfull_log, direction = "both", 
                        scope = list(lower = fitnull_log, upper = fitfull_log),
                        B = B, trans = "log",
                        y_col = "eqIncome", 
                        interval = interval,
                        trace = TRUE)

summary(fitopt_log$fit)
JcAIC_log <- fitopt_log$models[[length(fitopt_log$models)]]$JcAIC

opt_fixed_log <- fitopt_log$fit$call$fixed
#-------------------------------------------------------------------------------

# BoxCox Trafo------------------------------------------------------------------
summary(eusilcA_smp$eqIncome)
shift <- 0
# optimal lambda for the full model
lambdaopt <- optimal_parameter(generic_opt,
                              fixed = fixed,
                              smp_data = eusilcA_smp,
                              smp_domains = "district",
                              transformation = "box.cox",
                              interval = interval)
# Initial box-cox transformed y with the full model
eusilcA_smp$BC_eqIncome <- ((eusilcA_smp$eqIncome + shift)^lambdaopt - 1)/lambdaopt
fixed_BC <- update.formula(fixed, BC_eqIncome ~.)

fitfull_BC <- lme(fixed = fixed_BC, random = ~ 1|district, 
                  method = "ML", data = eusilcA_smp)
fitnull_BC <- lme(fixed = BC_eqIncome ~ 1, random = ~ 1|district, 
                  method = "ML", data = eusilcA_smp)

fitopt_BC <- stepJcAIC(object = fitfull_BC, direction = "both", 
                       scope = list(lower = fitnull_BC, upper = fitfull_BC),
                       B = B, trans = "box.cox",
                       y_col = "eqIncome", 
                       interval = interval,
                       trace = TRUE)

summary(fitopt_BC$fit)
JcAIC_BC <- fitopt_BC$models[[length(fitopt_BC$models)]]$JcAIC

opt_fixed_BC <- fitopt_BC$fit$call$fixed

#-------------------------------------------------------------------------------
# Optimal model
c(opt_fixed_no, opt_fixed_log, opt_fixed_BC)

which.min(c(JcAIC_no, JcAIC_log, JcAIC_BC))

################################################################################
#                             EBP with optimal model                           #
################################################################################
data("eusilcA_pop")
# B, L, MSE setting
L <- 100
B <- 50#0
MSE <- TRUE
# variables
opt_fixed_BC <- update.formula(opt_fixed_BC, eqIncome ~ .)

EBP_opt <- ebp(fixed = opt_fixed_BC, 
           pop_data = eusilcA_pop,
           pop_domains = "district", smp_data = eusilcA_smp,
           smp_domains = "district", transformation = "box.cox", 
           MSE = MSE, interval = c(-1, 1),
           na.rm = TRUE, boot_type = "wild", B = B, L = L)
summary(EBP_opt)
