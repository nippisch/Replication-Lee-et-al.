# analysis

load("data/test_r4_final_2.RData")
load("data/replication_R10_z.RData")
load("data/replication_R10_z_b50.RData")
library(dplyr)
library(stringr)

res_df <- data.frame(dgp = character(),
                     transformation = character(),
                     R = numeric(),
                     opt_model = character(),
                     x1 = numeric(),
                     x2 = numeric(),
                     x3 = numeric(),
                     #x4 = numeric(),
                     z1 = numeric(),
                     #z2 = numeric(),
                     opt_lambda = numeric(),
                     JcAIC = numeric(),
                     bc_boot = numeric())

# in case of empty list elements
sim_res_JcAIC <- lapply(sim_res_JcAIC, function(level2) {
  lapply(level2, function(level3) {
    Filter(Negate(is.null), level3)
  })
})

# res_norm_high <- sim_res_JcAIC[["normal_high"]]
# res_norm_low <- sim_res_JcAIC[["normal_low"]]
# res_log <- sim_res_JcAIC[["log"]]
# res_boxcox <- sim_res_JcAIC[["boxcox"]]
# 
# res_norm_high_no <- res_norm_high[["no"]]

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
      res_df[obs, "opt_model"] <- paste0(all.vars(res_tmp$opt_vars[[3]]), collapse = "+")
      res_df[obs, "x1"] <- ifelse(str_detect(res_df[obs, "opt_model"], "x1"), 1, 0)
      res_df[obs, "x2"] <- ifelse(str_detect(res_df[obs, "opt_model"], "x2"), 1, 0)
      res_df[obs, "x3"] <- ifelse(str_detect(res_df[obs, "opt_model"], "x3"), 1, 0)
      #res_df[obs, "x4"] <- ifelse(str_detect(res_df[obs, "opt_model"], "x4"), 1, 0)
      res_df[obs, "z1"] <- ifelse(str_detect(res_df[obs, "opt_model"], "z1"), 1, 0)
      #res_df[obs, "z2"] <- ifelse(str_detect(res_df[obs, "opt_model"], "z2"), 1, 0)
      res_df[obs, "opt_lambda"] <- res_tmp$lambdaopt
      res_df[obs, "JcAIC"] <- res_tmp$JcAIC
      
      obs <- obs + 1
      
    }
    
  }
  
}

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
  
