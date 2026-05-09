# =========================================================
# 2-Group (k=2) Blinded SSR Summaries & Figures
# Keep only:
#   1. sigma_tilde2       = naive blinded variance estimator
#   2. MSE                = unblinded variance estimator
#   3. adjusted_variance  = adjusted blinded variance estimator
#
# Add:
#   Table 4 = Bias of re-estimated per-group sample sizes
#   Table 5 = MSE  of re-estimated per-group sample sizes
# =========================================================

# ---- Settings ----
set.seed(20250826)

thetas        <- c(0.1, 0.3, 0.5, 0.7)
k_groups      <- 2
alpha         <- 0.05
target_power  <- 0.80
sigma_true    <- 10
mu_design     <- c(6, 1)   # planned means for 2 groups
out_dir       <- "outputs_ssr_g2"
R             <- 10000

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

theta_tag <- function(th) sprintf("theta%02d", as.integer(round(th * 100)))

# ---- Utility Functions ----

summarize_distribution <- function(x) {
  data.frame(
    Mean = mean(x, na.rm = TRUE),
    SD   = sd(x, na.rm = TRUE),
    P10  = as.numeric(quantile(x, 0.10, names = FALSE, na.rm = TRUE)),
    P25  = as.numeric(quantile(x, 0.25, names = FALSE, na.rm = TRUE)),
    P50  = as.numeric(quantile(x, 0.50, names = FALSE, na.rm = TRUE)),
    P75  = as.numeric(quantile(x, 0.75, names = FALSE, na.rm = TRUE)),
    P90  = as.numeric(quantile(x, 0.90, names = FALSE, na.rm = TRUE)),
    check.names = FALSE
  )
}

power_anova_equal_n <- function(n, k, alpha = 0.05, f2) {
  n <- pmax(2, round(n))
  df1 <- k - 1
  df2 <- k * (n - 1)
  lambda <- n * k * f2
  Fcrit <- qf(1 - alpha, df1 = df1, df2 = df2)
  pow <- 1 - pf(Fcrit, df1 = df1, df2 = df2, ncp = lambda)
  as.numeric(pow)
}

solve_n_exact <- function(k, f2, alpha = 0.05, target_power = 0.80,
                          n_min = 2, n_max = 50000) {
  if (!is.finite(f2) || f2 <= 0) return(NA_integer_)
  
  n <- max(2L, n_min)
  
  while (n <= n_max && power_anova_equal_n(n, k, alpha, f2) < target_power) {
    n <- ceiling(n * 1.2) + 1L
  }
  
  if (n > n_max) return(NA_integer_)
  
  low  <- max(2L, floor(n / 1.2))
  high <- n
  
  while (low < high) {
    mid <- floor((low + high) / 2)
    pwr <- power_anova_equal_n(mid, k, alpha, f2)
    
    if (pwr >= target_power) {
      high <- mid
    } else {
      low <- mid + 1L
    }
  }
  
  high
}

save_hist_pdf <- function(x_list, main_title, xlab, file_pdf) {
  pdf(file_pdf, width = 7, height = 5)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  }, add = TRUE)
  
  par(mfrow = c(length(x_list), 1), mar = c(4, 4, 3, 1))
  nm <- names(x_list)
  
  for (i in seq_along(x_list)) {
    hist(
      x_list[[i]],
      main = paste0(main_title, " (", nm[i], ")"),
      xlab = xlab,
      ylab = "Frequency"
    )
  }
}

save_boxplot_pdf <- function(x_list, main_title, ylab, file_pdf) {
  pdf(file_pdf, width = 7, height = 5)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  }, add = TRUE)
  
  boxplot(
    x_list,
    names = names(x_list),
    main = main_title,
    ylab = ylab
  )
}

# ---- Design helpers ----

effect_numerator <- function(mu_vec) {
  mu_bar <- mean(mu_vec)
  mean((mu_vec - mu_bar)^2)
}

initial_n_plan <- function(mu_vec, sigma, k, alpha = 0.05, target_power = 0.80) {
  num <- effect_numerator(mu_vec)
  f2  <- num / sigma^2
  solve_n_exact(k = k, f2 = f2, alpha = alpha, target_power = target_power, n_min = 2)
}

# ---- Interim variance estimators ----
# sigma_tilde2       = SSTO*/(N*-1)
# MSE                = SSE*/(N*-g)
# adjusted_variance  = [SSTO* - sum_i m_i*(mu_i-mu_w)^2]/(N*-g)

compute_interim_estimators <- function(y_mat, mu_vec, m_vec = NULL) {
  y_mat <- as.matrix(y_mat)
  k <- ncol(y_mat)
  
  if (length(mu_vec) != k) {
    stop("length(mu_vec) must equal ncol(y_mat).")
  }
  
  if (is.null(m_vec)) {
    m_vec <- rep(nrow(y_mat), k)
  }
  
  if (length(m_vec) != k) {
    stop("length(m_vec) must equal number of groups.")
  }
  
  N_star <- sum(m_vec)
  
  # collect data by group
  y_list <- lapply(seq_len(k), function(i) y_mat[seq_len(m_vec[i]), i])
  all_y <- unlist(y_list, use.names = FALSE)
  
  # SSTO*
  overall_mean <- mean(all_y)
  SSTO_star <- sum((all_y - overall_mean)^2)
  
  # naive blinded variance
  sigma_tilde2 <- SSTO_star / (N_star - 1)
  
  # unblinded variance
  group_means <- sapply(y_list, mean)
  SSE_star <- sum(sapply(seq_len(k), function(i) {
    sum((y_list[[i]] - group_means[i])^2)
  }))
  MSE <- SSE_star / (N_star - k)
  
  # adjusted blinded variance
  mu_bar_w <- sum((m_vec / N_star) * mu_vec)
  correction <- sum(m_vec * (mu_vec - mu_bar_w)^2)
  
  adjusted_variance <- (SSTO_star - correction) / (N_star - k)
  
  list(
    SSTO_star = SSTO_star,
    SSE_star = SSE_star,
    sigma_tilde2 = sigma_tilde2,
    MSE = MSE,
    adjusted_variance = adjusted_variance
  )
}

# ---- One simulation replicate ----
one_replicate_ssr <- function(theta, mu_design, sigma_true, alpha, target_power, k) {
  
  # Step 1: initial sample size (true/planned per-group sample size)
  n_plan <- initial_n_plan(
    mu_vec = mu_design,
    sigma = sigma_true,
    k = k,
    alpha = alpha,
    target_power = target_power
  )
  
  if (!is.finite(n_plan) || is.na(n_plan)) {
    return(list(
      n_plan = NA,
      sigma_tilde2 = NA,
      MSE = NA,
      adjusted_variance = NA,
      n_star = NA,
      n_star2 = NA,
      n_star3 = NA
    ))
  }
  
  # Step 2: interim sample size
  m <- max(2L, floor(theta * n_plan))
  m_vec <- rep(m, k)
  
  # Step 3: generate interim data
  y_mat <- sapply(mu_design, function(mu_i) {
    rnorm(m, mean = mu_i, sd = sigma_true)
  })
  
  if (!is.matrix(y_mat)) {
    y_mat <- matrix(y_mat, nrow = m, ncol = k)
  }
  
  # Step 4: compute estimators
  est <- compute_interim_estimators(
    y_mat = y_mat,
    mu_vec = mu_design,
    m_vec = m_vec
  )
  
  # Step 5: effect size re-estimation
  num <- effect_numerator(mu_design)
  
  f2_naive <- if (is.finite(est$sigma_tilde2) && est$sigma_tilde2 > 0) {
    num / est$sigma_tilde2
  } else {
    NA_real_
  }
  
  f2_mse <- if (is.finite(est$MSE) && est$MSE > 0) {
    num / est$MSE
  } else {
    NA_real_
  }
  
  f2_adj <- if (is.finite(est$adjusted_variance) && est$adjusted_variance > 0) {
    num / est$adjusted_variance
  } else {
    NA_real_
  }
  
  # Step 6: re-estimated per-group sample sizes
  n_star <- solve_n_exact(
    k = k,
    f2 = f2_naive,
    alpha = alpha,
    target_power = target_power,
    n_min = m
  )
  
  n_star2 <- solve_n_exact(
    k = k,
    f2 = f2_mse,
    alpha = alpha,
    target_power = target_power,
    n_min = m
  )
  
  n_star3 <- solve_n_exact(
    k = k,
    f2 = f2_adj,
    alpha = alpha,
    target_power = target_power,
    n_min = m
  )
  
  list(
    n_plan = n_plan,
    sigma_tilde2 = est$sigma_tilde2,
    MSE = est$MSE,
    adjusted_variance = est$adjusted_variance,
    n_star = n_star,
    n_star2 = n_star2,
    n_star3 = n_star3
  )
}

# ---- Main simulation storage ----
variance_results <- list()
n_results <- list()
power_results <- list()

# ---- Simulation loop ----
for (th in thetas) {
  tag <- theta_tag(th)
  
  sim_res <- replicate(R, {
    one_replicate_ssr(
      theta = th,
      mu_design = mu_design,
      sigma_true = sigma_true,
      alpha = alpha,
      target_power = target_power,
      k = k_groups
    )
  }, simplify = FALSE)
  
  n_plan_vec <- sapply(sim_res, function(x) x$n_plan)
  
  sigma_tilde2_vec <- sapply(sim_res, function(x) x$sigma_tilde2)
  MSE_vec <- sapply(sim_res, function(x) x$MSE)
  adjusted_variance_vec <- sapply(sim_res, function(x) x$adjusted_variance)
  
  n_star_vec  <- sapply(sim_res, function(x) x$n_star)
  n_star2_vec <- sapply(sim_res, function(x) x$n_star2)
  n_star3_vec <- sapply(sim_res, function(x) x$n_star3)
  
  variance_results[[tag]] <- list(
    sigma_tilde2 = sigma_tilde2_vec,
    MSE = MSE_vec,
    adjusted_variance = adjusted_variance_vec
  )
  
  n_results[[tag]] <- list(
    n_plan  = n_plan_vec,
    n_star  = n_star_vec,
    n_star2 = n_star2_vec,
    n_star3 = n_star3_vec
  )
  
  # true effect size for achieved power
  num <- effect_numerator(mu_design)
  f2_true <- num / sigma_true^2
  
  power_results[[tag]] <- list(
    n_star  = power_anova_equal_n(n_star_vec,  k = k_groups, alpha = alpha, f2 = f2_true),
    n_star2 = power_anova_equal_n(n_star2_vec, k = k_groups, alpha = alpha, f2 = f2_true),
    n_star3 = power_anova_equal_n(n_star3_vec, k = k_groups, alpha = alpha, f2 = f2_true)
  )
}

# =========================================================
# Table 2 and Table 3:
# Bias and MSE summaries for variance estimators
# =========================================================
for (th in thetas) {
  tag <- theta_tag(th)
  
  vres <- variance_results[[tag]]
  
  bias_mse_var_table <- data.frame(
    Estimator = c("sigma_tilde2", "MSE", "adjusted_variance"),
    Bias = c(
      mean(vres$sigma_tilde2, na.rm = TRUE) - sigma_true^2,
      mean(vres$MSE, na.rm = TRUE) - sigma_true^2,
      mean(vres$adjusted_variance, na.rm = TRUE) - sigma_true^2
    ),
    MSE = c(
      mean((vres$sigma_tilde2 - sigma_true^2)^2, na.rm = TRUE),
      mean((vres$MSE - sigma_true^2)^2, na.rm = TRUE),
      mean((vres$adjusted_variance - sigma_true^2)^2, na.rm = TRUE)
    )
  )
  
  write.csv(
    bias_mse_var_table,
    file.path(out_dir, paste0("table2_3_variance_bias_mse_", tag, ".csv")),
    row.names = FALSE
  )
}

# =========================================================
# Table 4:
# Bias of re-estimated sample size estimators
# =========================================================
for (th in thetas) {
  tag <- theta_tag(th)
  
  nres <- n_results[[tag]]
  
  valid_idx_1 <- is.finite(nres$n_plan) & is.finite(nres$n_star)
  valid_idx_2 <- is.finite(nres$n_plan) & is.finite(nres$n_star2)
  valid_idx_3 <- is.finite(nres$n_plan) & is.finite(nres$n_star3)
  
  table4_bias_n <- data.frame(
    Estimator = c("n_star", "n_star2", "n_star3"),
    Bias = c(
      mean(nres$n_star[valid_idx_1]  - nres$n_plan[valid_idx_1], na.rm = TRUE),
      mean(nres$n_star2[valid_idx_2] - nres$n_plan[valid_idx_2], na.rm = TRUE),
      mean(nres$n_star3[valid_idx_3] - nres$n_plan[valid_idx_3], na.rm = TRUE)
    )
  )
  
  write.csv(
    table4_bias_n,
    file.path(out_dir, paste0("table4_bias_n_", tag, ".csv")),
    row.names = FALSE
  )
}

# =========================================================
# Table 5:
# MSE of re-estimated sample size estimators
# =========================================================
for (th in thetas) {
  tag <- theta_tag(th)
  
  nres <- n_results[[tag]]
  
  valid_idx_1 <- is.finite(nres$n_plan) & is.finite(nres$n_star)
  valid_idx_2 <- is.finite(nres$n_plan) & is.finite(nres$n_star2)
  valid_idx_3 <- is.finite(nres$n_plan) & is.finite(nres$n_star3)
  
  table5_mse_n <- data.frame(
    Estimator = c("n_star", "n_star2", "n_star3"),
    MSE = c(
      mean((nres$n_star[valid_idx_1]  - nres$n_plan[valid_idx_1])^2, na.rm = TRUE),
      mean((nres$n_star2[valid_idx_2] - nres$n_plan[valid_idx_2])^2, na.rm = TRUE),
      mean((nres$n_star3[valid_idx_3] - nres$n_plan[valid_idx_3])^2, na.rm = TRUE)
    )
  )
  
  write.csv(
    table5_mse_n,
    file.path(out_dir, paste0("table5_mse_n_", tag, ".csv")),
    row.names = FALSE
  )
}

# =========================================================
# Table 6:
# Distribution summaries and plots
# =========================================================
for (th in thetas) {
  tag <- theta_tag(th)
  
  var_methods <- variance_results[[tag]]
  n_methods   <- n_results[[tag]]
  p_methods   <- power_results[[tag]]
  
  var_methods_clean <- lapply(var_methods, function(x) x[is.finite(x) & !is.na(x)])
  
  # exclude n_plan from Table 6 distribution summary
  n_methods_clean <- lapply(n_methods[c("n_star", "n_star2", "n_star3")], function(x) {
    x[is.finite(x) & !is.na(x)]
  })
  
  p_methods_clean <- lapply(p_methods, function(x) x[is.finite(x) & !is.na(x)])
  
  dist_var <- do.call(rbind, lapply(var_methods_clean, summarize_distribution))
  rownames(dist_var) <- names(var_methods_clean)
  
  dist_n <- do.call(rbind, lapply(n_methods_clean, summarize_distribution))
  rownames(dist_n) <- names(n_methods_clean)
  
  dist_power <- do.call(rbind, lapply(p_methods_clean, summarize_distribution))
  rownames(dist_power) <- names(p_methods_clean)
  
  write.csv(
    dist_var,
    file.path(out_dir, paste0("dist_variance_", tag, ".csv")),
    row.names = TRUE
  )
  write.csv(
    dist_n,
    file.path(out_dir, paste0("table6_dist_n_", tag, ".csv")),
    row.names = TRUE
  )
  write.csv(
    dist_power,
    file.path(out_dir, paste0("dist_power_", tag, ".csv")),
    row.names = TRUE
  )
  
  save_hist_pdf(
    x_list = var_methods_clean,
    main_title = paste0("Variance Estimators (", tag, ")"),
    xlab = "Estimated variance",
    file_pdf = file.path(out_dir, paste0("hist_variance_", tag, ".pdf"))
  )
  
  save_boxplot_pdf(
    x_list = var_methods_clean,
    main_title = paste0("Variance Estimators by Method (", tag, ")"),
    ylab = "Estimated variance",
    file_pdf = file.path(out_dir, paste0("box_variance_", tag, ".pdf"))
  )
  
  save_hist_pdf(
    x_list = n_methods_clean,
    main_title = paste0("Re-estimated Per-Group n (", tag, ")"),
    xlab = "Per-group sample size",
    file_pdf = file.path(out_dir, paste0("hist_n_", tag, ".pdf"))
  )
  
  save_boxplot_pdf(
    x_list = n_methods_clean,
    main_title = paste0("Per-Group n by Method (", tag, ")"),
    ylab = "Per-group sample size",
    file_pdf = file.path(out_dir, paste0("box_n_", tag, ".pdf"))
  )
  
  save_hist_pdf(
    x_list = p_methods_clean,
    main_title = paste0("Achieved Power (", tag, ")"),
    xlab = "Power",
    file_pdf = file.path(out_dir, paste0("hist_power_", tag, ".pdf"))
  )
}

cat("Done. Files written to: ", normalizePath(out_dir), "\n")

# =========================================================
# Print all main results to console
# =========================================================

cat("\n==============================\n")
cat("Main Simulation Results\n")
cat("==============================\n")

for (th in thetas) {
  tag <- theta_tag(th)
  
  cat("\n------------------------------------------\n")
  cat("Theta =", th, " (", tag, ")\n")
  cat("------------------------------------------\n")
  
  # ---- 1. Variance estimator Bias and MSE ----
  vres <- variance_results[[tag]]
  
  bias_mse_var_table <- data.frame(
    Estimator = c("sigma_tilde2", "MSE", "adjusted_variance"),
    Bias = c(
      mean(vres$sigma_tilde2, na.rm = TRUE) - sigma_true^2,
      mean(vres$MSE, na.rm = TRUE) - sigma_true^2,
      mean(vres$adjusted_variance, na.rm = TRUE) - sigma_true^2
    ),
    MSE = c(
      mean((vres$sigma_tilde2 - sigma_true^2)^2, na.rm = TRUE),
      mean((vres$MSE - sigma_true^2)^2, na.rm = TRUE),
      mean((vres$adjusted_variance - sigma_true^2)^2, na.rm = TRUE)
    )
  )
  
  cat("\nVariance Estimators: Bias and MSE\n")
  print(bias_mse_var_table, row.names = FALSE)
  
  # ---- 2. Distribution summaries for variance estimators ----
  var_methods_clean <- lapply(variance_results[[tag]], function(x) x[is.finite(x) & !is.na(x)])
  dist_var <- do.call(rbind, lapply(var_methods_clean, summarize_distribution))
  rownames(dist_var) <- names(var_methods_clean)
  
  cat("\nVariance Estimators: Distribution Summary\n")
  print(dist_var)
  
  # ---- 3. Bias of re-estimated sample sizes (Table 4) ----
  nres <- n_results[[tag]]
  
  valid_idx_1 <- is.finite(nres$n_plan) & is.finite(nres$n_star)
  valid_idx_2 <- is.finite(nres$n_plan) & is.finite(nres$n_star2)
  valid_idx_3 <- is.finite(nres$n_plan) & is.finite(nres$n_star3)
  
  table4_bias_n <- data.frame(
    Estimator = c("n_star", "n_star2", "n_star3"),
    Bias = c(
      mean(nres$n_star[valid_idx_1]  - nres$n_plan[valid_idx_1], na.rm = TRUE),
      mean(nres$n_star2[valid_idx_2] - nres$n_plan[valid_idx_2], na.rm = TRUE),
      mean(nres$n_star3[valid_idx_3] - nres$n_plan[valid_idx_3], na.rm = TRUE)
    )
  )
  
  cat("\nTable 4: Bias of Re-estimated Sample Sizes\n")
  print(table4_bias_n, row.names = FALSE)
  
  # ---- 4. MSE of re-estimated sample sizes (Table 5) ----
  table5_mse_n <- data.frame(
    Estimator = c("n_star", "n_star2", "n_star3"),
    MSE = c(
      mean((nres$n_star[valid_idx_1]  - nres$n_plan[valid_idx_1])^2, na.rm = TRUE),
      mean((nres$n_star2[valid_idx_2] - nres$n_plan[valid_idx_2])^2, na.rm = TRUE),
      mean((nres$n_star3[valid_idx_3] - nres$n_plan[valid_idx_3])^2, na.rm = TRUE)
    )
  )
  
  cat("\nTable 5: MSE of Re-estimated Sample Sizes\n")
  print(table5_mse_n, row.names = FALSE)
  
  # ---- 5. Distribution summaries for re-estimated sample sizes (Table 6) ----
  n_methods_clean <- lapply(n_results[[tag]][c("n_star", "n_star2", "n_star3")], function(x) {
    x[is.finite(x) & !is.na(x)]
  })
  dist_n <- do.call(rbind, lapply(n_methods_clean, summarize_distribution))
  rownames(dist_n) <- names(n_methods_clean)
  
  cat("\nTable 6: Re-estimated Per-Group Sample Sizes\n")
  print(dist_n)
  
  # ---- 6. Distribution summaries for achieved power ----
  p_methods_clean <- lapply(power_results[[tag]], function(x) x[is.finite(x) & !is.na(x)])
  dist_power <- do.call(rbind, lapply(p_methods_clean, summarize_distribution))
  rownames(dist_power) <- names(p_methods_clean)
  
  cat("\nAchieved Power Summary\n")
  print(dist_power)
}

cat("\n==============================\n")
cat("Done.\n")
cat("==============================\n")

# =========================================================
# Print key parameter values to console
# =========================================================
cat("\nKey Parameters:\n")
cat("thetas =", paste(thetas, collapse = ", "), "\n")
cat("k_groups =", k_groups, "\n")
cat("alpha =", alpha, "\n")
cat("target_power =", target_power, "\n")
cat("sigma_true =", sigma_true, "\n")
cat("mu1 =", mu_design[1], "\n")
cat("mu2 =", mu_design[2], "\n")

# =========================================================
# Optional: combine all Table 4 results into one CSV
# =========================================================
table4_all <- do.call(rbind, lapply(thetas, function(th) {
  tag <- theta_tag(th)
  nres <- n_results[[tag]]
  
  valid_idx_1 <- is.finite(nres$n_plan) & is.finite(nres$n_star)
  valid_idx_2 <- is.finite(nres$n_plan) & is.finite(nres$n_star2)
  valid_idx_3 <- is.finite(nres$n_plan) & is.finite(nres$n_star3)
  
  data.frame(
    Theta = th,
    Estimator = c("n_star", "n_star2", "n_star3"),
    Bias = c(
      mean(nres$n_star[valid_idx_1]  - nres$n_plan[valid_idx_1], na.rm = TRUE),
      mean(nres$n_star2[valid_idx_2] - nres$n_plan[valid_idx_2], na.rm = TRUE),
      mean(nres$n_star3[valid_idx_3] - nres$n_plan[valid_idx_3], na.rm = TRUE)
    )
  )
}))

write.csv(
  table4_all,
  file.path(out_dir, "table4_bias_n_all_theta.csv"),
  row.names = FALSE
)

# =========================================================
# Optional: combine all Table 5 results into one CSV
# =========================================================
table5_all <- do.call(rbind, lapply(thetas, function(th) {
  tag <- theta_tag(th)
  nres <- n_results[[tag]]
  
  valid_idx_1 <- is.finite(nres$n_plan) & is.finite(nres$n_star)
  valid_idx_2 <- is.finite(nres$n_plan) & is.finite(nres$n_star2)
  valid_idx_3 <- is.finite(nres$n_plan) & is.finite(nres$n_star3)
  
  data.frame(
    Theta = th,
    Estimator = c("n_star", "n_star2", "n_star3"),
    MSE = c(
      mean((nres$n_star[valid_idx_1]  - nres$n_plan[valid_idx_1])^2, na.rm = TRUE),
      mean((nres$n_star2[valid_idx_2] - nres$n_plan[valid_idx_2])^2, na.rm = TRUE),
      mean((nres$n_star3[valid_idx_3] - nres$n_plan[valid_idx_3])^2, na.rm = TRUE)
    )
  )
}))

write.csv(
  table5_all,
  file.path(out_dir, "table5_mse_n_all_theta.csv"),
  row.names = FALSE
)