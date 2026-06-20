############################################################
## Unequal Allocation Simulation for Exact ANOVA BSSR
############################################################

set.seed(12345)

## -----------------------------
## Exact ANOVA power and planning
## -----------------------------

anova_power_exact <- function(N, mu, sigma2, w, alpha = 0.05) {
  g  <- length(mu)
  df1 <- g - 1
  df2 <- N - g
  
  if (df2 <= 0) return(0)
  
  mu_w <- sum(w * mu)
  f2 <- sum(w * (mu - mu_w)^2) / sigma2
  lambda <- N * f2
  
  crit <- qf(1 - alpha, df1, df2)
  power <- 1 - pf(crit, df1, df2, ncp = lambda)
  
  return(power)
}

plan_N_exact <- function(mu, sigma2, w,
                         alpha = 0.05,
                         target_power = 0.80,
                         N_max = 5000) {
  g <- length(mu)
  
  for (N in seq(g + 2, N_max)) {
    if (anova_power_exact(N, mu, sigma2, w, alpha) >= target_power) {
      return(N)
    }
  }
  
  stop("N_max too small.")
}

## -----------------------------
## Allocate sample size by ratio
## -----------------------------

allocate_N <- function(N, w) {
  n <- floor(N * w)
  
  while (sum(n) < N) {
    k <- which.max(N * w - n)
    n[k] <- n[k] + 1
  }
  
  return(n)
}

## -----------------------------
## Variance estimators
## -----------------------------

variance_estimators <- function(y, group, mu_design, w_design) {
  N <- length(y)
  g <- length(mu_design)
  
  ybar <- mean(y)
  SSTO <- sum((y - ybar)^2)
  
  SSE <- sum(tapply(y, group, function(z) {
    sum((z - mean(z))^2)
  }))
  
  sigma2_NB <- SSTO / (N - 1)
  sigma2_UB <- SSE / (N - g)
  
  mu_w <- sum(w_design * mu_design)
  adj_term <- N * sum(w_design * (mu_design - mu_w)^2)
  
  sigma2_Adj <- (SSTO - adj_term) / (N - 1)
  sigma2_Adj <- max(sigma2_Adj, 1e-8)
  
  out <- c(
    NB  = sigma2_NB,
    UB  = sigma2_UB,
    Adj = sigma2_Adj
  )
  
  return(out)
}

## -----------------------------
## One simulation replicate
## -----------------------------

one_rep <- function(mu,
                    sigma,
                    w,
                    theta = 0.5,
                    alpha = 0.05,
                    target_power = 0.80) {
  
  g <- length(mu)
  sigma2 <- sigma^2
  
  N_star <- plan_N_exact(
    mu = mu,
    sigma2 = sigma2,
    w = w,
    alpha = alpha,
    target_power = target_power
  )
  
  N_int <- floor(theta * N_star)
  N_int <- max(N_int, g + 2)
  
  n_int <- allocate_N(N_int, w)
  
  group <- rep(seq_len(g), times = n_int)
  y <- unlist(lapply(seq_len(g), function(i) {
    rnorm(n_int[i], mean = mu[i], sd = sigma)
  }))
  
  sig2_hat <- variance_estimators(
    y = y,
    group = group,
    mu_design = mu,
    w_design = w
  )
  
  N_hat <- sapply(sig2_hat, function(s2) {
    plan_N_exact(
      mu = mu,
      sigma2 = s2,
      w = w,
      alpha = alpha,
      target_power = target_power
    )
  })
  
  N_final <- pmax(N_hat, N_int)
  
  achieved_power <- sapply(N_final, function(N) {
    anova_power_exact(
      N = N,
      mu = mu,
      sigma2 = sigma2,
      w = w,
      alpha = alpha
    )
  })
  
  data.frame(
    theta = theta,
    N_star = N_star,
    N_int = N_int,
    Method = names(sig2_hat),
    sigma2_hat = as.numeric(sig2_hat),
    N_hat = as.numeric(N_hat),
    N_final = as.numeric(N_final),
    Achieved_power = as.numeric(achieved_power)
  )
}

## -----------------------------
## Monte Carlo simulation
## -----------------------------

run_sim_unequal <- function(R = 5000,
                            mu = c(1, 4, 6),
                            sigma = 10,
                            w = c(0.25, 0.25, 0.50),
                            theta_vec = c(0.3, 0.5, 0.7),
                            alpha = 0.05,
                            target_power = 0.80) {
  
  results <- list()
  counter <- 1
  
  for (theta in theta_vec) {
    for (r in seq_len(R)) {
      results[[counter]] <- one_rep(
        mu = mu,
        sigma = sigma,
        w = w,
        theta = theta,
        alpha = alpha,
        target_power = target_power
      )
      counter <- counter + 1
    }
  }
  
  do.call(rbind, results)
}

## -----------------------------
## Run example
## -----------------------------

sim_unequal <- run_sim_unequal(
  R = 5000,
  mu = c(1, 4, 6),
  sigma = 10,
  w = c(0.25, 0.25, 0.50),
  theta_vec = c(0.3, 0.5, 0.7),
  alpha = 0.05,
  target_power = 0.80
)

## -----------------------------
## Summary table
## -----------------------------

summary_unequal <- aggregate(
  cbind(sigma2_hat, N_hat, N_final, Achieved_power) ~ theta + Method,
  data = sim_unequal,
  FUN = function(x) c(
    Mean = mean(x),
    SD = sd(x),
    Q25 = quantile(x, 0.25),
    Median = median(x),
    Q75 = quantile(x, 0.75)
  )
)

print(summary_unequal)

## Optional: cleaner summary
library(dplyr)

summary_table <- sim_unequal %>%
  group_by(theta, Method) %>%
  summarise(
    Mean_sigma2_hat = mean(sigma2_hat),
    SD_sigma2_hat   = sd(sigma2_hat),
    Mean_N_hat      = mean(N_hat),
    SD_N_hat        = sd(N_hat),
    Median_N_hat    = median(N_hat),
    Q25_N_hat       = quantile(N_hat, 0.25),
    Q75_N_hat       = quantile(N_hat, 0.75),
    Mean_power      = mean(Achieved_power),
    SD_power        = sd(Achieved_power),
    .groups = "drop"
  )

print(summary_table)
