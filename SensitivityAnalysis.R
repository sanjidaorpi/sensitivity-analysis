future::plan("multisession", workers = 5)



# Set the file paths
here::i_am("sensitivity_analysis/Testset.R")

# Load the functions for the causal testbed and generate some data
# Generating a prime number just so that we can check against the recycle rule
source(here::here("sensitivity_analysis", "CausalTestbed.R"))

# Load the genetic operators
source(here::here("sensitivity_analysis", "GeneticOperators.R"))



# ----- Get initial estimate -----



# We use the data generating process to generate observational data for 10000 people
set.seed(344029)
sample_obs <- obs_dgp(n = 500)

# Create an IPW model
vaccine_model <- glm(
  vaccine ~ sex + age,
  data = sample_obs,
  family = "binomial"
)

# Estimate propensity
prob_V1 <- predict(vaccine_model, newdata = sample_obs, type = "response")
prob_V0 <- (1 - prob_V1)
prob_vaccine <- ifelse(sample_obs$vaccine == 1, prob_V1, prob_V0)

# We measure the effect of getting and not getting a vaccine using IPW
V_0 <- ifelse(sample_obs$vaccine == 0, 1, 0)
CE_0 <- sum((1 / prob_vaccine) * sample_obs$outcome * V_0) / sum(V_0)
V_1 <- ifelse(sample_obs$vaccine == 1, 1, 0)
CE_1 <- sum((1 / prob_vaccine) * sample_obs$outcome * V_1) / sum(V_1)

# So that the IPW corrected average treatment effect is
ate_ipw <- (CE_1 - CE_0)
ate_ipw



# ----- Bootstrap Sensitivity Analysis -----



# Start the sensitivity analysis
G = 1.1
ate_maximized = FALSE

while (!ate_maximized) {
  
  # Compute bootstrap estimates of ATE
  ate_list <- future.apply::future_lapply(1:100, function(i){
    boot_inds <- sample(1:nrow(sample_obs), nrow(sample_obs), replace = TRUE)
    boot_data <- sample_obs[boot_inds,]
    
    vaccine_model <- glm(
      vaccine ~ sex + age,
      data = boot_data,
      family = "binomial"
    )
    
    prob_V1 <- predict(vaccine_model, newdata = boot_data, type = "response")
    prob_V0 <- (1 - prob_V1)
    prob_vaccine <- ifelse(boot_data$vaccine == 1, prob_V1, prob_V0)
    
    a_val <- genetic_optimizer(boot_data, prob_vaccine, G)
    V_0 <- ifelse(boot_data$vaccine == 0, 1, 0)
    CE_0 <- sum((a_val / prob_vaccine) * boot_data$outcome * V_0) / sum(V_0)
    V_1 <- ifelse(boot_data$vaccine == 1, 1, 0)
    CE_1 <- sum((a_val / prob_vaccine) * boot_data$outcome * V_1) / sum(V_1)
    boot_ipw <- (CE_1 - CE_0)
    
    return(boot_ipw)
  }, future.seed = TRUE)
  
  # Bootstrap Confidence Interval
  ate_list <- unlist(ate_list)
  delta_stars <- (ate_list - ate_ipw)
  d = quantile(delta_stars, c(0.025, 0.975))
  ci = ate_ipw - c(d[2], d[1])
  
  # Check if zero is inside the confidence interval
  if ((ci[1] < 0) & (ci[2] > 0)) {
    break
  } else {
    G <- G + 0.1
  }
  
}

# Naive Bayes
nb.fit <- e1071::naiveBayes(vaccine ~ sex + age, data = sample_obs, type = "raw")
response_nb <- predict(nb.fit, newdata = sample_obs, type = "raw")
prob_V1_nb <- response_nb[ ,2]
#table(sample_obs$vaccine, response_nb)
prob_V0 <- (1 - prob_V1_nb)
prob_vaccine <- ifelse(sample_obs$vaccine == 1, prob_V1_nb, prob_V0)
V_0 <- ifelse(sample_obs$vaccine == 0, 1, 0)
CE_0 <- sum((1 / prob_vaccine) * sample_obs$outcome * V_0) / sum(V_0)
V_1 <- ifelse(sample_obs$vaccine == 1, 1, 0)
CE_1 <- sum((1 / prob_vaccine) * sample_obs$outcome * V_1) / sum(V_1)
ate_ipw <- (CE_1 - CE_0)

# Linear Regression
prob_V1_lr <- predict(vaccine_model, newdata = sample_obs, type = "response")
response_lr <- ifelse(prob_V1_lr >= 0.5, 1, 0)
#table(sample_obs$vaccine, response_lr)
prob_V0 <- (1 - prob_V1_lr)
prob_vaccine <- ifelse(sample_obs$vaccine == 1, prob_V1_lr, prob_V0)
V_0 <- ifelse(sample_obs$vaccine == 0, 1, 0)
CE_0 <- sum((1 / prob_vaccine) * sample_obs$outcome * V_0) / sum(V_0)
V_1 <- ifelse(sample_obs$vaccine == 1, 1, 0)
CE_1 <- sum((1 / prob_vaccine) * sample_obs$outcome * V_1) / sum(V_1)
ate_ipw <- (CE_1 - CE_0)

future::plan("sequential")
