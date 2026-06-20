############################################################
## Full Revised R Code
## Unequal Allocation + Non-normal Robustness Simulation
## Exact ANOVA Planning and Adjusted Blinded Variance Estimation
############################################################

set.seed(12345)

############################################################
## 1. Exact ANOVA power
############################################################

anova_power_exact <- function(N,
                              mu,
                              sigma2,
                              w,
                              alpha = 0.05) {
  
  g <- length(mu)
  df1 <- g - 1
  df2 <- N - g
  
  if (df2 <= 0) return(0)
  if (is.na(sigma2) || sigma2 <= 0) return(NA)
  
  mu_w <- sum(w * mu)
  f2 <- sum(w * (mu - mu_w)^2) / sigma2
  lambda <- N * f2
  
  crit <- qf(1 - alpha, df1, df2)
  
  power <- 1 - pf(
    crit,
    df1 = df1,
    df2 = df2,
    ncp = lambda
  )
  
  return(power)
}

############################################################
## 2. Exact sample size planning
## Revised: larger N_max and no stopping error
############################################################

plan_N_exact <- function(mu,
                         sigma2,
                         w,
                         alpha = 0.05,
                         target_power = 0.80,
                         N_max = 100000) {
  
  g <- length(mu)
  
  if (is.na(sigma2) || sigma2 <= 0) {
    return(NA)
  }
  
  for (N in seq(g + 2, N_max)) {
    
    pow <- anova_power_exact(
      N = N,
      mu = mu,
      sigma2 = sigma2,
      w = w,
      alpha = alpha
    )
    
    if (!is.na(pow) && pow >= target_power) {
      return(N)
    }
  }
  
  return(NA)
}

############################################################
## 3. Group sample sizes under unequal allocation
############################################################

get_group_sizes <- function(N, w) {
  
  g <- length(w)
  n_i <- floor(N * w)
  n_i[g] <- N - sum(n_i[-g])
  
  n_i <- pmax(n_i, 2)
  
  diff_N <- N - sum(n_i)
  
  if (diff_N != 0) {
    n_i[g] <- n_i[g] + diff_N
  }
  
  return(n_i)
}

############################################################
## 4. Generate data under normal and non-normal distributions
############################################################

generate_data_nonnormal <- function(N,
                                    mu,
                                    sigma,
                                    w,
                                    dist = "normal") {
  
  g <- length(mu)
  n_i <- get_group_sizes(N, w)
  
  y <- numeric()
  group <- numeric()
  
  for (i in 1:g) {
    
    ni <- n_i[i]
    
    if (dist == "normal") {
      
      eps <- rnorm(ni, mean = 0, sd = sigma)
      
    } else if (dist == "t3") {
      
      eps <- sigma * rt(ni, df = 3) / sqrt(3)
      
    } else if (dist == "t5") {
      
      eps <- sigma * rt(ni, df = 5) / sqrt(5 / 3)
      
    } else if (dist == "lognormal") {
      
      raw <- rlnorm(ni, meanlog = 0, sdlog = 0.6)
      raw <- raw - mean(raw)
      raw <- raw / sd(raw)
      eps <- sigma * raw
      
    } else if (dist == "gamma") {
      
      raw <- rgamma(ni, shape = 2, scale = 1)
      raw <- raw - mean(raw)
      raw <- raw / sd(raw)
      eps <- sigma * raw
      
    } else if (dist == "contaminated") {
      
      contam <- rbinom(ni, size = 1, prob = 0.05)
      eps <- rnorm(
        ni,
        mean = 0,
        sd = sigma * ifelse(contam == 1, 5, 1)
      )
      
    } else {
      
      stop("Unknown distribution.")
    }
    
    y <- c(y, mu[i] + eps)
    group <- c(group, rep(i, ni))
  }
  
  data.frame(
    y = y,
    group = factor(group)
  )
}

############################################################
## 5. Variance estimators
## NB  = naive blinded
## UB  = unblinded
## Adj = adjusted blinded
############################################################

variance_estimators <- function(dat,
                                mu,
                                w) {
  
  y <- dat$y
  group <- dat$group
  
  N <- length(y)
  g <- length(mu)
  
  SSTO <- sum((y - mean(y))^2)
  
  SSE <- sum(
    tapply(y, group, function(x) {
      sum((x - mean(x))^2)
    })
  )
  
  n_i <- as.numeric(table(group))
  w_i <- n_i / sum(n_i)
  mu_w <- sum(w_i * mu)
  
  planned_between <- sum(n_i * (mu - mu_w)^2)
  
  sigma2_NB <- SSTO / (N - 1)
  sigma2_UB <- SSE / (N - g)
  sigma2_Adj <- (SSTO - planned_between) / (N - 1)
  
  sigma2_Adj <- max(sigma2_Adj, 1e-8)
  
  c(
    NB = sigma2_NB,
    UB = sigma2_UB,
    Adj = sigma2_Adj
  )
}

############################################################
## 6. Final ANOVA test
############################################################

anova_reject <- function(dat,
                         alpha = 0.05) {
  
  fit <- aov(y ~ group, data = dat)
  pval <- summary(fit)[[1]][["Pr(>F)"]][1]
  
  as.numeric(pval < alpha)
}

############################################################
## 7. Empirical achieved power
############################################################

estimate_power_empirical <- function(N,
                                     mu,
                                     sigma,
                                     w,
                                     dist,
                                     alpha = 0.05,
                                     B_power = 300) {
  
  reject <- numeric(B_power)
  
  for (b in 1:B_power) {
    
    dat <- generate_data_nonnormal(
      N = N,
      mu = mu,
      sigma = sigma,
      w = w,
      dist = dist
    )
    
    reject[b] <- anova_reject(dat, alpha = alpha)
  }
  
  mean(reject)
}

############################################################
## 8. Main simulation function
############################################################

run_sim_unequal_nonnormal <- function(R = 1000,
                                      mu = c(10, 13, 15),
                                      sigma = 10,
                                      w = c(0.25, 0.25, 0.50),
                                      theta_vec = c(0.3, 0.5, 0.7),
                                      dist_vec = c(
                                        "normal",
                                        "t3",
                                        "t5",
                                        "lognormal",
                                        "gamma",
                                        "contaminated"
                                      ),
                                      alpha = 0.05,
                                      target_power = 0.80,
                                      B_power = 300,
                                      N_max = 100000) {
  
  sigma2 <- sigma^2
  g <- length(mu)
  
  N_star <- plan_N_exact(
    mu = mu,
    sigma2 = sigma2,
    w = w,
    alpha = alpha,
    target_power = target_power,
    N_max = N_max
  )
  
  if (is.na(N_star)) {
    stop("Initial N_star could not be found. Increase N_max or use larger effect size.")
  }
  
  all_results <- data.frame()
  
  for (dist in dist_vec) {
    
    for (theta in theta_vec) {
      
      N_interim <- max(g + 2, floor(theta * N_star))
      
      for (estimator in c("NB", "UB", "Adj")) {
        
        cat(
          "Running:",
          "dist =", dist,
          "| theta =", theta,
          "| estimator =", estimator,
          "\n"
        )
        
        sigma2_hat_vec <- numeric(R)
        N_hat_vec <- numeric(R)
        N_final_vec <- numeric(R)
        power_vec <- numeric(R)
        failed_vec <- numeric(R)
        
        for (r in 1:R) {
          
          dat_int <- generate_data_nonnormal(
            N = N_interim,
            mu = mu,
            sigma = sigma,
            w = w,
            dist = dist
          )
          
          s2_all <- variance_estimators(
            dat = dat_int,
            mu = mu,
            w = w
          )
          
          s2 <- s2_all[estimator]
          
          N_hat <- plan_N_exact(
            mu = mu,
            sigma2 = s2,
            w = w,
            alpha = alpha,
            target_power = target_power,
            N_max = N_max
          )
          
          if (is.na(N_hat)) {
            failed_vec[r] <- 1
            N_hat <- N_max
          }
          
          N_final <- max(N_hat, N_interim)
          
          emp_power <- estimate_power_empirical(
            N = N_final,
            mu = mu,
            sigma = sigma,
            w = w,
            dist = dist,
            alpha = alpha,
            B_power = B_power
          )
          
          sigma2_hat_vec[r] <- s2
          N_hat_vec[r] <- N_hat
          N_final_vec[r] <- N_final
          power_vec[r] <- emp_power
        }
        
        tmp <- data.frame(
          Distribution = dist,
          Theta = theta,
          Estimator = estimator,
          N_star = N_star,
          N_interim = N_interim,
          Mean_sigma2_hat = mean(sigma2_hat_vec, na.rm = TRUE),
          Bias_sigma2_hat = mean(sigma2_hat_vec, na.rm = TRUE) - sigma2,
          MSE_sigma2_hat = mean((sigma2_hat_vec - sigma2)^2, na.rm = TRUE),
          Mean_N_hat = mean(N_hat_vec, na.rm = TRUE),
          SD_N_hat = sd(N_hat_vec, na.rm = TRUE),
          Median_N_hat = median(N_hat_vec, na.rm = TRUE),
          Q25_N_hat = quantile(N_hat_vec, 0.25, na.rm = TRUE),
          Q75_N_hat = quantile(N_hat_vec, 0.75, na.rm = TRUE),
          Mean_N_final = mean(N_final_vec, na.rm = TRUE),
          Achieved_Power = mean(power_vec, na.rm = TRUE),
          Power_Deviation = mean(power_vec, na.rm = TRUE) - target_power,
          Failure_Rate = mean(failed_vec)
        )
        
        all_results <- rbind(all_results, tmp)
      }
    }
  }
  
  all_results$Mean_sigma2_hat <- round(all_results$Mean_sigma2_hat, 3)
  all_results$Bias_sigma2_hat <- round(all_results$Bias_sigma2_hat, 3)
  all_results$MSE_sigma2_hat <- round(all_results$MSE_sigma2_hat, 3)
  
  all_results$Mean_N_hat <- round(all_results$Mean_N_hat, 2)
  all_results$SD_N_hat <- round(all_results$SD_N_hat, 2)
  all_results$Median_N_hat <- round(all_results$Median_N_hat, 2)
  all_results$Q25_N_hat <- round(all_results$Q25_N_hat, 2)
  all_results$Q75_N_hat <- round(all_results$Q75_N_hat, 2)
  all_results$Mean_N_final <- round(all_results$Mean_N_final, 2)
  
  all_results$Achieved_Power <- round(all_results$Achieved_Power, 3)
  all_results$Power_Deviation <- round(all_results$Power_Deviation, 3)
  all_results$Failure_Rate <- round(all_results$Failure_Rate, 3)
  
  return(all_results)
}

############################################################
## 9. Run simulation
############################################################

sim_nonnormal <- run_sim_unequal_nonnormal(
  R = 1000,
  mu = c(10, 13, 15),
  sigma = 10,
  w = c(0.25, 0.25, 0.50),
  theta_vec = c(0.3, 0.5, 0.7),
  dist_vec = c(
    "normal",
    "t3",
    "t5",
    "lognormal",
    "gamma",
    "contaminated"
  ),
  alpha = 0.05,
  target_power = 0.80,
  B_power = 300,
  N_max = 100000
)

print(sim_nonnormal)

############################################################
## 10. Export results
############################################################

write.csv(
  sim_nonnormal,
  "bssr_unequal_nonnormal_robustness_results.csv",
  row.names = FALSE
)
