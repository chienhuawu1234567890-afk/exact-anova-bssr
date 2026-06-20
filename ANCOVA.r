############################################################
## ANCOVA Simulation for Exact ANOVA BSSR
## R_Z^2 treated as a design-stage parameter
############################################################

set.seed(12345)

## -----------------------------
## Exact ANCOVA power
## -----------------------------

ancova_power_exact <- function(N, mu, sigma2, w,
                               R2z = 0.30,
                               alpha = 0.05) {
  g <- length(mu)

  df1 <- g - 1
  df2 <- N - g - 1   # one baseline covariate

  if (df2 <= 0) return(0)

  mu_w <- sum(w * mu)

  sigma2_ancova <- (1 - R2z) * sigma2

  f2_ancova <- sum(w * (mu - mu_w)^2) / sigma2_ancova

  lambda <- N * f2_ancova

  crit <- qf(1 - alpha, df1, df2)

  power <- 1 - pf(crit, df1, df2, ncp = lambda)

  return(power)
}

## -----------------------------
## Exact ANCOVA sample size planning
## -----------------------------

plan_N_ancova_exact <- function(mu, sigma2, w,
                                R2z = 0.30,
                                alpha = 0.05,
                                target_power = 0.80,
                                N_max = 5000) {
  g <- length(mu)

  for (N in seq(g + 2, N_max)) {
    power_N <- ancova_power_exact(
      N = N,
      mu = mu,
      sigma2 = sigma2,
      w = w,
      R2z = R2z,
      alpha = alpha
    )

    if (power_N >= target_power) {
      return(N)
    }
  }

  return(NA)
}

## -----------------------------
## Blinded variance estimators
## -----------------------------

variance_estimators_blinded <- function(Y, mu, w, theta) {
  N_I <- length(Y)

  mu_w <- sum(w * mu)

  delta2 <- sum(w * (mu - mu_w)^2)

  ## Naive blinded variance
  sigma2_NB <- var(Y)

  ## Design-adjusted blinded variance
  ## theta = N_I / N, so N_I / N = theta
  sigma2_Adj <- sigma2_NB - theta * delta2

  sigma2_Adj <- max(sigma2_Adj, 1e-8)

  return(list(
    NB = sigma2_NB,
    Adj = sigma2_Adj
  ))
}

## -----------------------------
## Unblinded ANCOVA residual variance estimator
## Used only as a benchmark
## -----------------------------

variance_estimator_unblinded_ancova <- function(Y, A, Z) {
  fit <- lm(Y ~ factor(A) + Z)

  sigma2_UB <- sum(resid(fit)^2) / df.residual(fit)

  return(sigma2_UB)
}

## -----------------------------
## One simulation replicate
## -----------------------------

simulate_one_ancova <- function(mu,
                                sigma2 = 100,
                                w = c(1/3, 1/3, 1/3),
                                theta = 0.50,
                                R2z = 0.30,
                                alpha = 0.05,
                                target_power = 0.80) {
  g <- length(mu)

  ## Planned total sample size using true sigma2
  N_plan <- plan_N_ancova_exact(
    mu = mu,
    sigma2 = sigma2,
    w = w,
    R2z = R2z,
    alpha = alpha,
    target_power = target_power
  )

  N_I <- ceiling(theta * N_plan)

  n_I <- as.vector(rmultinom(1, size = N_I, prob = w))

  A <- rep(seq_len(g), times = n_I)

  ## Generate baseline covariate
  Z <- rnorm(N_I, mean = 0, sd = 1)

  ## Choose gamma so that R_Z^2 is approximately the design value
  ## R_Z^2 = gamma^2 Var(Z) / {gamma^2 Var(Z) + sigma_e^2}
  ## Here Var(Z)=1 and residual error variance is sigma2.
  gamma <- sqrt(R2z * sigma2 / (1 - R2z))

  eps <- rnorm(N_I, mean = 0, sd = sqrt(sigma2))

  Y <- mu[A] + gamma * Z + eps

  ## Blinded variance estimators on total pooled outcome
  est_blind <- variance_estimators_blinded(
    Y = Y,
    mu = mu,
    w = w,
    theta = theta
  )

  ## Convert blinded outcome variance estimators to ANCOVA residual variance
  sigma2_NB_ANCOVA  <- (1 - R2z) * est_blind$NB
  sigma2_Adj_ANCOVA <- (1 - R2z) * est_blind$Adj

  ## Unblinded ANCOVA benchmark
  sigma2_UB_ANCOVA <- variance_estimator_unblinded_ancova(
    Y = Y,
    A = A,
    Z = Z
  )

  ## Re-estimated sample sizes
  N_NB <- plan_N_ancova_exact(
    mu = mu,
    sigma2 = sigma2_NB_ANCOVA,
    w = w,
    R2z = R2z,
    alpha = alpha,
    target_power = target_power
  )

  N_Adj <- plan_N_ancova_exact(
    mu = mu,
    sigma2 = sigma2_Adj_ANCOVA,
    w = w,
    R2z = R2z,
    alpha = alpha,
    target_power = target_power
  )

  N_UB <- plan_N_ancova_exact(
    mu = mu,
    sigma2 = sigma2_UB_ANCOVA,
    w = w,
    R2z = R2z,
    alpha = alpha,
    target_power = target_power
  )

  ## Final sample size cannot be smaller than interim sample size
  N_final_NB  <- max(N_I, N_NB)
  N_final_Adj <- max(N_I, N_Adj)
  N_final_UB  <- max(N_I, N_UB)

  ## Achieved power under true ANCOVA model
  power_NB <- ancova_power_exact(
    N = N_final_NB,
    mu = mu,
    sigma2 = sigma2,
    w = w,
    R2z = R2z,
    alpha = alpha
  )

  power_Adj <- ancova_power_exact(
    N = N_final_Adj,
    mu = mu,
    sigma2 = sigma2,
    w = w,
    R2z = R2z,
    alpha = alpha
  )

  power_UB <- ancova_power_exact(
    N = N_final_UB,
    mu = mu,
    sigma2 = sigma2,
    w = w,
    R2z = R2z,
    alpha = alpha
  )

  out <- data.frame(
    theta = theta,
    Method = c("NB", "Adj", "UB"),
    R2z = R2z,
    N_plan = N_plan,
    N_interim = N_I,

    sigma2_hat = c(
      sigma2_NB_ANCOVA,
      sigma2_Adj_ANCOVA,
      sigma2_UB_ANCOVA
    ),

    N_hat = c(
      N_NB,
      N_Adj,
      N_UB
    ),

    N_final = c(
      N_final_NB,
      N_final_Adj,
      N_final_UB
    ),

    Achieved_power = c(
      power_NB,
      power_Adj,
      power_UB
    )
  )

  return(out)
}

## -----------------------------
## Run simulation
## -----------------------------

run_sim_ancova <- function(R = 5000,
                           mu = c(1, 4, 6),
                           sigma = 10,
                           w = c(1/3, 1/3, 1/3),
                           theta_vec = c(0.3, 0.5, 0.7),
                           R2z_vec = c(0.10, 0.30, 0.50),
                           alpha = 0.05,
                           target_power = 0.80) {
  res <- list()
  counter <- 1

  sigma2 <- sigma^2

  for (R2z in R2z_vec) {
    for (theta in theta_vec) {
      for (r in seq_len(R)) {
        res[[counter]] <- simulate_one_ancova(
          mu = mu,
          sigma2 = sigma2,
          w = w,
          theta = theta,
          R2z = R2z,
          alpha = alpha,
          target_power = target_power
        )

        counter <- counter + 1
      }
    }
  }

  sim <- do.call(rbind, res)

  return(sim)
}

## -----------------------------
## Example
## -----------------------------

sim_ancova <- run_sim_ancova(
  R = 5000,
  mu = c(1, 4, 6),
  sigma = 10,
  w = c(1/3, 1/3, 1/3),
  theta_vec = c(0.3, 0.5, 0.7),
  R2z_vec = c(0.10, 0.30, 0.50),
  alpha = 0.05,
  target_power = 0.80
)

## -----------------------------
## Summary table
## -----------------------------

summary_ancova <- aggregate(
  cbind(sigma2_hat, N_hat, N_final, Achieved_power) ~
    R2z + theta + Method,
  data = sim_ancova,
  FUN = function(x) c(
    Mean = mean(x),
    SD = sd(x),
    Q25 = quantile(x, 0.25),
    Median = median(x),
    Q75 = quantile(x, 0.75)
  )
)

summary_ancova <- do.call(data.frame, summary_ancova)

names(summary_ancova) <- gsub("\\.", "_", names(summary_ancova))

print(summary_ancova)

## -----------------------------
## Save results
## -----------------------------

write.csv(
  sim_ancova,
  file = "ancova_bssr_simulation_raw.csv",
  row.names = FALSE
)

write.csv(
  summary_ancova,
  file = "ancova_bssr_simulation_summary.csv",
  row.names = FALSE
)
