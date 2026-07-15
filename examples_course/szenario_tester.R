#########################################################################
# Advanced Small Area Estimation 
#
# Author:    Prof. Dr. Timo Schmid
# Content:   Szenario tester              
#########################################################################

library(emdi)
library(ggplot2)
library(knitr)

# 1. Compare the density of the target variable in the population with 
# that of a single Monte Carlo sample.

test_smp <- Pop[[1]]

# Extract the population data from the sample
test_pop <- attr(test_smp, "pop")

# Calculate the threshold, mean, and median
test_thold <- 0.6 * median(test_smp$y)
test_mean <- mean(test_smp$y)
test_median <- median(test_smp$y)


ggplot() +
  geom_density(data = test_smp, aes(x = y, fill = "Sample"), alpha = 0.5) +
  geom_density(data = test_pop, aes(x = y, fill = "Population"), alpha = 0.5) +
  geom_vline(xintercept = test_thold, linetype = "dashed", color = "red", size = 1) +
  geom_vline(xintercept = test_mean, linetype = "dashed", color = "green", size = 1) +
  geom_vline(xintercept = test_median, linetype = "dashed", color = "blue", size = 1) +
  labs(title = "",
       x = "Target variable",
       y = "Density",
       fill = "Group") +
  scale_fill_manual(values = c("Sample" = "skyblue", "Population" = "orange")) +
  theme_minimal() +
  theme(legend.position = "top")


# 2. Estimate an EBP model with box cox transformation and look at a summary.
ebpBC <- ebp(y~x, pop_data = test_pop, pop_domains = "idD",
             smp_data = test_smp, smp_domains = "idD", L = 50, 
             transformation = "box.cox")
summary(ebpBC)


# 3. Plot the EBP model and look at the residuals. 
plot(ebpBC)


# 4. Take a closer look at some point estimates and check whether 
# the results make sense. 
tableSum = rbind(summary(ebpBC$ind$Mean),
                 summary(ebpBC$ind$Head_Count),
                 summary(ebpBC$ind$Poverty_Gap),
                 summary(ebpBC$ind$Quintile_Share))
rownames(tableSum) = c("Mean","HCR", "PGAP", "QSR")
colnames(tableSum) = c("Min", "Q1", "Median", "Mean", "Q3", "Max")
kable(tableSum, digits=4, caption ="Summary of indicators")


# 5. Calculate the transformation parameter for the entire populations
# and examine the distribution.
all_ebps <- vector("list", length(Pop))
for(i in seq_along(Pop))
{
  tmp_smp <- Pop[[i]]
  tmp_pop <- attr(tmp_smp, "pop")  
  tmp_thold <- 0.6 * median(tmp_smp$y)
  
  all_ebps[[i]] <-ebp(y~x, pop_data = tmp_pop, pop_domains = "idD",
                      smp_data = tmp_smp, smp_domains = "idD", L = 2, 
                      transformation = "box.cox")
  
}
# Extract the optimal_lambda from each ebp result
all_lambdas <- sapply(all_ebps, function(list) list[["transform_param"]][["optimal_lambda"]])

# Combine the lambdas into a data frame
lambda_data <- data.frame(population = rep(1:length(all_ebps), each = length(all_lambdas[[1]])),
                          optimal_lambda = unlist(all_lambdas))

ggplot(lambda_data, aes(y = optimal_lambda)) +
  geom_boxplot(fill = "skyblue", color = "black") +
  labs(title = "Optimal transformation parameter across populations",
       x = "",
       y = "Transformation parameter") +
  theme_minimal()

