
library(bcf)
library(grf)
library(MASS)
library(mvtnorm)
library(dplyr)
library(doParallel)
library(foreach)

set.seed(2020)
n_sim <- 50
n <- 250
p <- 5
results_file <- "bcf_simulation_results.csv"

if (file.exists(results_file)) file.remove(results_file)

#g(x4) mapping for categorical variable
g_fun <- function(x4) {
  vals <- c(2, -1, -4)
  return(vals[x4])
}

#treatment effect function
tau_fun <- list(
  homogeneous = function(x) rep(3, nrow(x)),
  heterogeneous = function(x) 1 + 2 * x[,2] * x[,5]
)

#prognostic function
mu_fun <- list(
  linear = function(x) 1 + g_fun(x[,4]) + x[,1] * x[,3],
  nonlinear = function(x) -6 + g_fun(x[,4]) + 6 * abs(x[,3] - 1)
)

#propensity function
propensity_score <- function(mu, x1) {
  s <- sd(mu)
  prob <- 0.8 * pnorm(3 * mu/s - 0.5 * x1) + 0.05 + runif(length(mu)) / 10
  return(pmin(pmax(prob, 0.01), 0.99))
}

#estimation metrics
compute_metrics <- function(estimates, true_values, ci_lower, ci_upper) {
  rmse <- sqrt(mean((estimates - true_values)^2))
  coverage <- mean(ci_lower <= true_values & true_values <= ci_upper)
  len <- mean(ci_upper - ci_lower)
  return(c(rmse=rmse, coverage=coverage, len=len))
}

#PARALLELIZATION. USE IN COMPUTE CLUSTER WITH MULTIPLE CORES IF POSSIBLE. 
#note: 48 cores and 256gb memory takes ~11 hours in 'Oscar' to run

#detects available cores (adjust as needed for your cluster setup)
n_cores <- min(parallel::detectCores(), 48)  # Use up to 48 cores
cat(sprintf("Setting up parallel processing with %d cores\n", n_cores))

#creates and registers cluster - IMPORTANT: must be done before starting the loop
cl <- makeCluster(n_cores)

#exports necessary functions and variables to all worker nodes
clusterExport(cl, c("g_fun", "tau_fun", "mu_fun", "propensity_score", "compute_metrics"))

#loads required packages on all nodes
clusterEvalQ(cl, {
  library(bcf)
  library(MASS)
  library(mvtnorm)
  library(dplyr)
})

#registers the parallel backend
registerDoParallel(cl)

cat("Starting parallel simulations...\n")

#creates a dataframe with all simulation configurations
sim_grid <- expand.grid(
  sim = 1:n_sim,
  mu_type = names(mu_fun),
  tau_type = names(tau_fun),
  stringsAsFactors = FALSE
)


results <- foreach(i = 1:nrow(sim_grid), .combine = rbind, 
                   .packages = c("bcf", "MASS", "mvtnorm", "dplyr")) %dopar% {
                     
     current_sim <- sim_grid[i, "sim"]
     current_mu_type <- as.character(sim_grid[i, "mu_type"])
     current_tau_type <- as.character(sim_grid[i, "tau_type"])
     
     cat(sprintf("Running simulation %d: mu=%s, tau=%s\n", 
                 current_sim, current_mu_type, current_tau_type))

      #generate X
      x1 <- rnorm(n)
      x2 <- rnorm(n)
      x3 <- rnorm(n)
      x4 <- sample(1:3, n, replace = TRUE)
      x5 <- rbinom(n, 1, 0.5)
      X <- cbind(x1, x2, x3, x4, x5)
      
      x4_oh <- model.matrix(~ as.factor(x4) - 1)
      X_model <- cbind(x1, x2, x3, x5, x4_oh)
      
      #accesses functions through explicit namespace to avoid scope issues
      #gets the correct function from mu_fun and tau_fun
      mu_function <- mu_fun[[current_mu_type]]
      tau_function <- tau_fun[[current_tau_type]]
      
      mu <- mu_function(X)
      tau <- tau_function(X)
      pi <- propensity_score(mu, x1)
      z <- rbinom(n, 1, pi)
      y <- mu + tau * z + rnorm(n)
      
      fit_bcf <- bcf(
        y = y,
        z = z,
        x_control = X_model,
        x_moderate = X_model,
        pihat = pi,
        nburn = 1000,
        nsim = 2000,
        verbose = FALSE
      )
      
      tau_post <- fit_bcf$tau
      est_ate <- rowMeans(tau_post)
      est_ate_mean <- mean(est_ate)
      est_ate_ci <- quantile(est_ate, probs = c(0.025, 0.975))
      
      true_ate <- mean(tau)
      ate_metrics <- compute_metrics(est_ate_mean, true_ate, est_ate_ci[1], est_ate_ci[2])
      
      est_cate_mean <- colMeans(tau_post)
      est_cate_ci_lwr <- apply(tau_post, 2, quantile, 0.025)
      est_cate_ci_upr <- apply(tau_post, 2, quantile, 0.975)
      cate_metrics <- compute_metrics(est_cate_mean, tau, est_cate_ci_lwr, est_cate_ci_upr)
      
      data.frame(
        method = "bcf",
        dgp_mu = current_mu_type,
        dgp_tau = current_tau_type,
        sim = current_sim,
        ate_rmse = ate_metrics["rmse"],
        ate_cover = ate_metrics["coverage"],
        ate_len = ate_metrics["len"],
        cate_rmse = cate_metrics["rmse"],
        cate_cover = cate_metrics["coverage"],
        cate_len = cate_metrics["len"]
      )
 }
      
      
write.csv(results, file = results_file, row.names = FALSE)

#stops the cluster when done
stopCluster(cl)
cat("All simulations completed!\n")

#ggplot trace plot
plot = ggplot(tau_long, aes(x = Iteration, y = Tau, color = Unit)) +
  geom_line() +
  labs(title = "Trace Plots for 3 Random Units (BCF)",
       x = "MCMC Iteration",
       y = "Posterior Draws of Tau") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")

      
  