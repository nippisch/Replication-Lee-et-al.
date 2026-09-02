# analysis

load("data/test_r4_final_2.RData")
load("data/replication_R10_z.RData")
load("data/replication_R10_z_b50.RData")
load("data/simulation_study_results_raw.RData")
load("data/replication_study_results_raw.RData")
load("/Volumes/HP v212w/simulation_study_results_MSE_raw.RData")
library(dplyr)
library(stringr)

res_df <- data.frame(dgp = character(),
                     transformation = character(),
                     R = numeric(),
                     opt_model = character(),
                     x1 = numeric(),
                     x2 = numeric(),
                     x3 = numeric(),
                     z1 = numeric(),
                     z2 = numeric(),
                     opt_lambda = numeric(),
                     JcAIC = numeric(),
                     boot_BC = numeric(),
                     clogLik = numeric(),
                     logJacobian = numeric(),
                     MSE = numeric())

obs <- 1

for (i in 1:length(sim_res_JcAIC)) {
  
  res_dgp <- sim_res_JcAIC[[i]]
  
  for (j in 1:length(res_dgp)) {
    
    res_trafo <- res_dgp[[j]]
    
    for (r in 1:length(res_trafo)) {
      
      res_tmp <- res_trafo[[r]]
      
      res_df[obs, "dgp"] <- names(sim_res_JcAIC)[i]
      res_df[obs, "transformation"] <- names(res_dgp)[j]
      res_df[obs, "R"] <- r
      if(names(res_tmp[1]) == "error") {
        res_df[obs, c("opt_model", "x1", "x2", "x3", "z1", "z2", "opt_lambda", "JcAIC")] <- NA
        obs <- obs + 1
        next
      }
      res_df[obs, "opt_model"] <- paste0(all.vars(res_tmp$opt_vars[[3]]), collapse = "+")
      res_df[obs, "x1"] <- ifelse(str_detect(res_df[obs, "opt_model"], "x1"), 1, 0)
      res_df[obs, "x2"] <- ifelse(str_detect(res_df[obs, "opt_model"], "x2"), 1, 0)
      res_df[obs, "x3"] <- ifelse(str_detect(res_df[obs, "opt_model"], "x3"), 1, 0)
      #res_df[obs, "x4"] <- ifelse(str_detect(res_df[obs, "opt_model"], "x4"), 1, 0)
      res_df[obs, "z1"] <- ifelse(str_detect(res_df[obs, "opt_model"], "z1"), 1, 0)
      #res_df[obs, "z2"] <- ifelse(str_detect(res_df[obs, "opt_model"], "z2"), 1, 0)
      res_df[obs, "opt_lambda"] <- ifelse(is.null(res_tmp$lambdaopt), res_tmp$JcAIC_result[["lambdaopt"]], res_tmp$lambdaopt)
      res_df[obs, "JcAIC"] <- res_tmp$JcAIC
      res_df[obs, "MSE"] <- res_tmp$MSE
      if(!is.null(res_tmp$JcAIC_result)) {
        res_df[obs, "boot_BC"] <- res_tmp$JcAIC_result$bootBC
        res_df[obs, "clogLik"] <- res_tmp$JcAIC_result$clogLik
        res_df[obs, "logJacobian"] <- res_tmp$JcAIC_result$logJacobian
      }
      
      
      obs <- obs + 1
      
    }
    
  }
  
}

res_df_replication <- res_df |> 
  select(-z2)

save(res_df_replication, file = "data/replication_study_results_processed.RData")

save(res_df, file = "data/simulation_study_results_processed.RData")

save(res_df, file = "data/simulation_study_results_MSE_processed.RData")

res_df_b50 <- res_df 
res_df_b50 |> 
  mutate(correct_model = if_else(x1 == 1 & x2 == 1 & x3 == 1 & z1 == 0, 1, 0)) |> 
  group_by(dgp, transformation) |> 
  summarise(mean_JcAIC = mean(JcAIC, na.rm = TRUE), x1 = sum(x1), x2 = sum(x2), x3 = sum(x3), z1 = sum(z1), mean_lambda = mean(opt_lambda, na.rm = TRUE),
            proportion_correct = mean(correct_model, na.rm = TRUE))

save(res_df, file = "data/results_replication_R10_z.RData")

res_df |> 
  mutate(correct_model = if_else(x1 == 1 & x2 == 1 & x3 == 1 & z1 == 0 & z2 == 0, 1, 0)) |> 
  group_by(dgp, transformation) |> 
  summarise(mean_JcAIC = mean(JcAIC, na.rm = TRUE), x1 = sum(x1), x2 = sum(x2), x3 = sum(x3), z1 = sum(z1), z2 = sum(z2), mean_lambda = mean(opt_lambda, na.rm = TRUE),
            proportion_correct = mean(correct_model, na.rm = TRUE))

load("data/results_replication_R10_z.RData")

res_df |> 
  mutate(correct_model = if_else(x1 == 1 & x2 == 1 & x3 == 1 & z1 == 0, 1, 0)) |> 
  group_by(dgp, transformation) |> 
  summarise(mean_JcAIC = mean(JcAIC, na.rm = TRUE), x1 = sum(x1), x2 = sum(x2), x3 = sum(x3), z1 = sum(z1), mean_lambda = mean(opt_lambda, na.rm = TRUE),
            proportion_correct = mean(correct_model, na.rm = TRUE))
  

load("data/sim_lambda_1_75.RData")
sim_res_JcAIC1 <- sim_res_JcAIC
sim_res_sample1 <- sim_res_sample
load("data/sim_lambda_151_225.RData")
sim_res_JcAIC3 <- sim_res_JcAIC
sim_res_sample3 <- sim_res_sample
load("data/sim_lambda_226_300.RData")
sim_res_JcAIC4 <- sim_res_JcAIC
sim_res_sample4 <- sim_res_sample

rm(sim_res_sample, sim_res_JcAIC)

sim_res_JcAIC <- list(
  normal_low = list(
    box.cox = do.call(
      c,
      lapply(
        list(sim_res_JcAIC1, sim_res_JcAIC2, sim_res_JcAIC3, sim_res_JcAIC4),
        function(x) x$normal_low$box.cox
      )
    )
  ),
  normal_high = list(
    box.cox = do.call(
      c,
      lapply(
        list(sim_res_JcAIC1, sim_res_JcAIC2, sim_res_JcAIC3, sim_res_JcAIC4),
        function(x) x$normal_high$box.cox
      )
    )
  )
)



res_df_lambda <- data.frame(dgp = character(),
                     transformation = character(),
                     R = numeric(),
                     opt_model = character(),
                     x1 = numeric(),
                     x2 = numeric(),
                     x3 = numeric(),
                     #x4 = numeric(),
                     z1 = numeric(),
                     z2 = numeric(),
                     opt_lambda = numeric(),
                     JcAIC = numeric(),
                     boot_BC = numeric(),
                     clogLik = numeric(),
                     lambda_start = numeric())

obs <- 1

for (i in 1:length(sim_res_JcAIC)) {
  
  res_dgp <- sim_res_JcAIC[[i]]
  
  for (j in 1:length(res_dgp)) {
    
    res_trafo <- res_dgp[[j]]
    
    for (r in 1:length(res_trafo)) {
      
      res_tmp <- res_trafo[[r]]
      
      res_df_lambda[obs, "dgp"] <- names(sim_res_JcAIC)[i]
      res_df_lambda[obs, "transformation"] <- names(res_dgp)[j]
      res_df_lambda[obs, "R"] <- r
      if(names(res_tmp[1]) == "error") {
        res_df_lambda[obs, c("opt_model", "x1", "x2", "x3", "z1", "z2", "opt_lambda", "JcAIC")] <- NA
        obs <- obs + 1
        next
      }
      res_df_lambda[obs, "opt_model"] <- paste0(all.vars(res_tmp$opt_vars[[3]]), collapse = "+")
      res_df_lambda[obs, "x1"] <- ifelse(str_detect(res_df_lambda[obs, "opt_model"], "x1"), 1, 0)
      res_df_lambda[obs, "x2"] <- ifelse(str_detect(res_df_lambda[obs, "opt_model"], "x2"), 1, 0)
      res_df_lambda[obs, "x3"] <- ifelse(str_detect(res_df_lambda[obs, "opt_model"], "x3"), 1, 0)
      res_df_lambda[obs, "z1"] <- ifelse(str_detect(res_df_lambda[obs, "opt_model"], "z1"), 1, 0)
      res_df_lambda[obs, "z2"] <- ifelse(str_detect(res_df_lambda[obs, "opt_model"], "z2"), 1, 0)
      res_df_lambda[obs, "opt_lambda"] <- res_tmp$lambdaopt
      res_df_lambda[obs, "JcAIC"] <- res_tmp$JcAIC
      if(!is.null(res_tmp$JcAIC_result)) {
        res_df_lambda[obs, "boot_BC"] <- res_tmp$JcAIC_result$bootBC
        res_df_lambda[obs, "clogLik"] <- res_tmp$JcAIC_result$clogLik
        res_df_lambda[obs, "lambda_start"] <- res_tmp$lambdaopt_start
      }
      
      
      obs <- obs + 1
      
    }
    
  }
  
}

res_df_lambda$difflambda <- res_df_lambda$lambda_start - res_df_lambda$opt_lambda
save(res_df_lambda, file = "data/sim_lambda_res_all.RData")

dens_df <- map_dfr(
  seq_along(sim_res_sample$normal_high),
  function(i) {
    d <- density(sim_res_sample$normal_high[[i]]$y)
    
    data.frame(
      x = d$x,
      density = d$y,
      sample = factor(i)
    )
  }
)
save(dens_df, file = "data/density_df.RData")
