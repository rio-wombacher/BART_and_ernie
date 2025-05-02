
library(bcf)
library(grf)
library(MASS)
library(mvtnorm)
library(dplyr)

set.seed(2020)
n_sim <- 50
n <- 250
p <- 5
results_file <- "bcf_grf_simulation_results.csv"

if (file.exists(results_file)) file.remove(results_file)

# g(x4) mapping for categorical variable
g_fun <- function(x4) {
  vals <- c(2, -1, -4)
  return(vals[x4])
}

# Treatment effect function
tau_fun <- list(
  homogeneous = function(x) rep(3, nrow(x)),
  heterogeneous = function(x) 1 + 2 * x[,2] * x[,5]
)

# Prognostic function
mu_fun <- list(
  linear = function(x) 1 + g_fun(x[,4]) + x[,1] * x[,3],
  nonlinear = function(x) -6 + g_fun(x[,4]) + 6 * abs(x[,3] - 1)
)

# Propensity function
propensity_score <- function(mu, x1) {
  s <- sd(mu)
  prob <- 0.8 * pnorm(3 * mu / s - 0.5 * x1) + 0.05 + runif(length(mu)) / 10
  return(pmin(pmax(prob, 0.01), 0.99))
}

# Estimation metrics
compute_metrics <- function(estimates, true_values, ci_lower, ci_upper) {
  rmse <- sqrt(mean((estimates - true_values)^2))
  coverage <- mean(ci_lower <= true_values & true_values <= ci_upper)
  len <- mean(ci_upper - ci_lower)
  return(c(rmse=rmse, coverage=coverage, len=len))
}

# Progress bar
total_steps <- length(mu_fun) * length(tau_fun) * n_sim * 2  # 2 methods
pb <- txtProgressBar(min = 0, max = total_steps, style = 3)
step <- 0

# Main loop
for (sim in 1:n_sim) {
  for (mu_type in names(mu_fun)) {
    for (tau_type in names(tau_fun)) {
      # Generate X
      x1 <- rnorm(n)
      x2 <- rnorm(n)
      x3 <- rnorm(n)
      x4 <- sample(1:3, n, replace = TRUE)
      x5 <- rbinom(n, 1, 0.5)
      X <- cbind(x1, x2, x3, x4, x5)
      
      x4_oh <- model.matrix(~ as.factor(x4) - 1)
      X_model <- cbind(x1, x2, x3, x5, x4_oh)
      
      mu <- mu_fun[[mu_type]](X)
      tau <- tau_fun[[tau_type]](X)
      pi <- propensity_score(mu, x1)
      z <- rbinom(n, 1, pi)
      y <- mu + tau * z + rnorm(n)
      
      ## --- Method 1: BCF ---
      fit_bcf <- bcf(
        y = y,
        z = z,
        x_control = X_model,
        x_moderate = X_model,
        pihat = pi,
        nburn = 200,
        nsim = 200,
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
      
      out_bcf <- data.frame(
        method = "bcf",
        dgp_mu = mu_type,
        dgp_tau = tau_type,
        sim = sim,
        ate_rmse = ate_metrics["rmse"],
        ate_cover = ate_metrics["coverage"],
        ate_len = ate_metrics["len"],
        cate_rmse = cate_metrics["rmse"],
        cate_cover = cate_metrics["coverage"],
        cate_len = cate_metrics["len"]
      )
      
      write.table(out_bcf, file = results_file, sep = ",", row.names = FALSE,
                  col.names = !file.exists(results_file), append = TRUE)
      step <- step + 1
      setTxtProgressBar(pb, step)
      
      # ## --- Method 2: Causal Forest ---
      # cf <- causal_forest(X_model, y, z)
      # cf_preds <- predict(cf, estimate.variance = TRUE)
      # 
      # est_cate_cf <- cf_preds$predictions
      # se_cate_cf <- sqrt(cf_preds$variance.estimates)
      # ci_lwr <- est_cate_cf - 1.96 * se_cate_cf
      # ci_upr <- est_cate_cf + 1.96 * se_cate_cf
      # 
      # cate_metrics_cf <- compute_metrics(est_cate_cf, tau, ci_lwr, ci_upr)
      # 
      # est_ate_cf <- mean(est_cate_cf)
      # se_ate_cf <- sd(est_cate_cf) / sqrt(n)
      # ate_ci_cf <- c(est_ate_cf - 1.96 * se_ate_cf, est_ate_cf + 1.96 * se_ate_cf)
      # ate_metrics_cf <- compute_metrics(est_ate_cf, true_ate, ate_ci_cf[1], ate_ci_cf[2])
      # 
      # out_cf <- data.frame(
      #   method = "grf",
      #   dgp_mu = mu_type,
      #   dgp_tau = tau_type,
      #   sim = sim,
      #   ate_rmse = ate_metrics_cf["rmse"],
      #   ate_cover = ate_metrics_cf["coverage"],
      #   ate_len = ate_metrics_cf["len"],
      #   cate_rmse = cate_metrics_cf["rmse"],
      #   cate_cover = cate_metrics_cf["coverage"],
      #   cate_len = cate_metrics_cf["len"]
      # )
      # 
      # write.table(out_cf, file = results_file, sep = ",", row.names = FALSE,
      #             col.names = FALSE, append = TRUE)
      # step <- step + 1
      # setTxtProgressBar(pb, step)
    }
  }
}
close(pb)
